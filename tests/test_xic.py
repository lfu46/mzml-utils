"""Tests for the xic module."""

import numpy as np
import pytest

from mzml_utils import Spectrum
from mzml_utils.constants import NEUTRON_MASS
from mzml_utils.xic import (extract_xic, extract_xics, XIC, ChromatogramSet,
                            coelution_score, isotope_trace_evidence)


def _ms1(scan, rt, peaks):
    """peaks: list of (mz, intensity)."""
    mz = np.array([p[0] for p in peaks], dtype=float)
    inten = np.array([p[1] for p in peaks], dtype=float)
    return Spectrum(scan_num=scan, mz=mz, intensity=inten, ms_level=1, rt=rt,
                    tic=float(inten.sum()),
                    base_peak_intensity=float(inten.max()) if len(inten) else 0.0)


def _ms2(scan, rt, peaks, precursor_mz=0.0, activation="HCD"):
    mz = np.array([p[0] for p in peaks], dtype=float)
    inten = np.array([p[1] for p in peaks], dtype=float)
    return Spectrum(scan_num=scan, mz=mz, intensity=inten, ms_level=2, rt=rt,
                    precursor_mz=precursor_mz, activation_type=activation,
                    filter_string=f"FTMS + p NSI d Full ms2 {precursor_mz}@{activation.lower()}",
                    tic=float(inten.sum()))


class FakeReader:
    """Minimal reader exposing iter_spectra(), like MzMLReader / SpectrumCache."""

    def __init__(self, specs):
        self._specs = specs

    def iter_spectra(self):
        for s in self._specs:
            yield s


@pytest.fixture
def ms1_reader():
    # Target ion 1000.50 present in 3 MS1 scans (apex at rt=11), a decoy peak, and an MS2 scan.
    specs = [
        _ms1(1, 10.0, [(500.0, 50.0), (1000.50, 100.0)]),
        _ms1(3, 11.0, [(1000.50, 500.0), (1200.0, 20.0)]),
        _ms1(5, 12.0, [(1000.50, 200.0)]),
        _ms2(2, 10.5, [(204.087, 8000.0), (138.055, 3000.0)], precursor_mz=1000.50),
    ]
    return FakeReader(specs)


class TestExtractXic:
    def test_basic_trace(self, ms1_reader):
        xic = extract_xic(ms1_reader, 1000.50, ms_level=1, tolerance=20.0, unit="ppm")
        assert isinstance(xic, XIC)
        # Only the 3 MS1 scans contribute; sorted by RT.
        assert xic.n_points == 3
        np.testing.assert_allclose(xic.rt, [10.0, 11.0, 12.0])
        np.testing.assert_allclose(xic.intensity, [100.0, 500.0, 200.0])

    def test_apex_and_nl(self, ms1_reader):
        xic = extract_xic(ms1_reader, 1000.50, ms_level=1)
        assert xic.max_intensity == pytest.approx(500.0)   # the "NL" value
        assert xic.apex_rt == pytest.approx(11.0)
        assert xic.area > 0

    def test_ms_level_filter_excludes_ms2(self, ms1_reader):
        # 204.087 lives only in the MS2 scan; an ms_level=1 XIC must not see it.
        xic = extract_xic(ms1_reader, 204.087, ms_level=1)
        assert xic.max_intensity == 0.0

    def test_ms2_fragment_xic(self, ms1_reader):
        xic = extract_xic(ms1_reader, 204.087, ms_level=2)
        assert xic.n_points == 1
        assert xic.intensity[0] == pytest.approx(8000.0)

    def test_rt_range(self, ms1_reader):
        xic = extract_xic(ms1_reader, 1000.50, ms_level=1, rt_range=(10.5, 12.5))
        assert xic.n_points == 2
        np.testing.assert_allclose(xic.rt, [11.0, 12.0])

    def test_tolerance_too_tight_misses(self, ms1_reader):
        # target off by ~0.5 Da -> 1 ppm window misses everything
        xic = extract_xic(ms1_reader, 1001.0, ms_level=1, tolerance=1.0, unit="ppm")
        assert xic.max_intensity == 0.0

    def test_da_unit(self, ms1_reader):
        xic = extract_xic(ms1_reader, 1000.50, ms_level=1, tolerance=0.02, unit="Da")
        assert xic.n_points == 3

    def test_auto_name(self, ms1_reader):
        xic = extract_xic(ms1_reader, 1000.50, ms_level=1)
        assert xic.name == "1000.5000"


class TestExtractXics:
    def test_multi_target_shared_axis(self, ms1_reader):
        cs = extract_xics(ms1_reader, {"decoy": 500.0, "target": 1000.50}, ms_level=1)
        assert isinstance(cs, ChromatogramSet)
        assert set(cs.names()) == {"decoy", "target"}
        # shared RT axis across lanes
        np.testing.assert_allclose(cs["target"].rt, cs["decoy"].rt)
        assert cs.n_scans == 3
        # 500.0 only in first MS1 scan
        np.testing.assert_allclose(cs["decoy"].intensity, [50.0, 0.0, 0.0])

    def test_tic_and_base_peak_collected(self, ms1_reader):
        cs = extract_xics(ms1_reader, [1000.50], ms_level=1, collect_tic=True)
        assert cs.tic.shape == cs.rt.shape
        assert cs.base_peak.shape == cs.rt.shape
        assert cs.tic[1] == pytest.approx(520.0)   # 500 + 20 at rt=11

    def test_mode_sum_vs_max(self):
        # two peaks within a wide 0.1 Da window -> sum=300, max=200
        spec = _ms1(1, 5.0, [(700.00, 100.0), (700.05, 200.0)])
        reader = FakeReader([spec])
        s = extract_xic(reader, 700.02, ms_level=1, tolerance=0.1, unit="Da", mode="sum")
        m = extract_xic(reader, 700.02, ms_level=1, tolerance=0.1, unit="Da", mode="max")
        assert s.intensity[0] == pytest.approx(300.0)
        assert m.intensity[0] == pytest.approx(200.0)

    def test_precursor_filter(self):
        specs = [
            _ms2(1, 5.0, [(204.087, 100.0)], precursor_mz=800.0),
            _ms2(2, 6.0, [(204.087, 900.0)], precursor_mz=1000.5),
        ]
        reader = FakeReader(specs)
        cs = extract_xics(reader, [204.087], ms_level=2,
                          precursor_mz=1000.5, precursor_tol=0.7)
        assert cs.n_scans == 1
        assert cs["204.0870"].intensity[0] == pytest.approx(900.0)

    def test_activation_filter(self):
        specs = [
            _ms2(1, 5.0, [(204.087, 100.0)], activation="HCD"),
            _ms2(2, 6.0, [(204.087, 900.0)], activation="ETD"),
        ]
        reader = FakeReader(specs)
        cs = extract_xics(reader, [204.087], ms_level=2, activation="HCD")
        assert cs.n_scans == 1
        assert cs["204.0870"].intensity[0] == pytest.approx(100.0)

    def test_empty_reader(self):
        cs = extract_xics(FakeReader([]), [204.087], ms_level=1)
        assert cs.n_scans == 0
        assert cs["204.0870"].max_intensity == 0.0


class TestPeakShape:
    """fwhm / points_across_peak: chromatographic width from the same object that
    already carries apex_rt and area, so no project re-derives it.

    The metric they serve is the instrument's own "desired minimum points across the
    peak" setting -- measured-versus-spec is what says whether the chromatography is
    adequately sampled, and whether a longer gradient has room to broaden peaks.
    """

    @staticmethod
    def _xic(intensity, rt=None):
        inten = np.asarray(intensity, dtype=float)
        rt = np.arange(len(inten), dtype=float) if rt is None else np.asarray(rt, dtype=float)
        return XIC(target_mz=500.0, rt=rt, intensity=inten)

    def test_fwhm_interpolates_between_scans(self):
        # half-max = 5; crossings interpolate to rt 2.25 and 5.75.
        # Snapping to the nearest scan instead would quantise the width in steps of a
        # whole cycle time -- the same order as the differences being measured.
        x = self._xic([0, 2, 4, 8, 10, 8, 4, 2, 0])
        assert x.fwhm == pytest.approx(3.5)
        assert x.points_across_peak == 3
        assert not x.fwhm_is_truncated
        assert x.apex_rt == pytest.approx(4.0)

    def test_second_peak_does_not_widen_the_first(self):
        # An XIC routinely holds a second elution or a co-eluting interference at the
        # same m/z; counting every point above half max anywhere would merge them.
        x = self._xic([0, 10, 0, 0, 0, 0, 9, 0])
        assert x.points_across_peak == 1
        assert x.fwhm < 1.5

    def test_truncated_peak_is_flagged_not_silently_reported(self):
        assert self._xic([2, 4, 8, 10]).fwhm_is_truncated       # apex at the last point
        assert self._xic([10, 8, 4, 2]).fwhm_is_truncated       # apex at the first
        assert not self._xic([0, 10, 0]).fwhm_is_truncated

    def test_degenerate_traces_return_zero_not_an_exception(self):
        for trace in ([], [0.0, 0.0, 0.0], [5.0]):
            x = self._xic(trace)
            assert x.fwhm == 0.0
            assert x.points_across_peak in (0, 1)

    def test_fwhm_uses_real_retention_times(self):
        narrow = self._xic([0, 10, 0], rt=[10.0, 10.2, 10.4])
        wide = self._xic([0, 10, 0], rt=[10.0, 12.0, 14.0])
        assert 0.0 < narrow.fwhm < 0.4
        assert wide.fwhm > narrow.fwhm

    def test_flat_top_peak_spans_the_plateau(self):
        x = self._xic([0, 6, 10, 10, 10, 6, 0])
        assert x.points_across_peak == 5       # three 10s plus both 6s (>= half max)
        assert x.fwhm > 3.0


class FilteringReader(FakeReader):
    """Reader whose iter_spectra takes the push-down filters, like the built-in readers."""

    def __init__(self, specs):
        super().__init__(specs)
        self.calls = []

    def iter_spectra(self, ms_level=None, rt_range=None):
        self.calls.append(dict(ms_level=ms_level, rt_range=rt_range))
        for s in self._specs:
            if ms_level is not None and s.ms_level != ms_level:
                continue
            if rt_range is not None and not (rt_range[0] <= s.rt <= rt_range[1]):
                continue
            yield s


class TestFilterPushDown:
    def test_filters_are_passed_to_a_reader_that_accepts_them(self, ms1_reader):
        reader = FilteringReader(ms1_reader._specs)
        extract_xic(reader, 1000.50, ms_level=1, rt_range=(10.5, 12.5))
        assert reader.calls == [dict(ms_level=1, rt_range=(10.5, 12.5))]

    def test_result_is_identical_with_and_without_push_down(self, ms1_reader):
        plain = extract_xic(ms1_reader, 1000.50, ms_level=1, rt_range=(10.5, 12.5))
        pushed = extract_xic(FilteringReader(ms1_reader._specs), 1000.50,
                             ms_level=1, rt_range=(10.5, 12.5))
        np.testing.assert_allclose(pushed.rt, plain.rt)
        np.testing.assert_allclose(pushed.intensity, plain.intensity)

    def test_none_filters_are_not_pushed(self, ms1_reader):
        reader = FilteringReader(ms1_reader._specs)
        extract_xics(reader, [1000.50], ms_level=None)
        assert reader.calls == [dict(ms_level=None, rt_range=None)]


def _trace(values, name="t"):
    rt = np.arange(len(values), dtype=float)
    return XIC(target_mz=500.0, rt=rt, intensity=np.asarray(values, dtype=float), name=name)


class TestCoelutionScore:
    def test_proportional_traces_score_one(self):
        a = _trace([0, 10, 50, 100, 40, 5, 0])
        assert coelution_score(a, _trace(a.intensity * 0.3)) == pytest.approx(1.0)

    def test_disjoint_traces_score_zero(self):
        assert coelution_score(_trace([0, 100, 50, 0, 0, 0]),
                               _trace([0, 0, 0, 0, 80, 40])) == pytest.approx(0.0)

    def test_partial_overlap_is_a_true_cosine(self):
        a, b = _trace([0, 3, 4, 0]), _trace([0, 0, 4, 3])
        assert coelution_score(a, b) == pytest.approx(16.0 / 25.0)

    def test_rt_range_restricts_the_comparison(self):
        a, b = _trace([100, 0, 10, 20, 10]), _trace([0, 100, 5, 10, 5])
        assert coelution_score(a, b) < 0.1
        assert coelution_score(a, b, rt_range=(2.0, 4.0)) == pytest.approx(1.0)

    def test_empty_or_flat_trace_scores_zero(self):
        assert coelution_score(_trace([0, 0, 0]), _trace([1, 2, 3])) == 0.0

    def test_mismatched_axes_raise(self):
        with pytest.raises(ValueError):
            coelution_score(_trace([1, 2, 3]), _trace([1, 2]))


def _envelope_reader(mono_mz, charge, ratios, extra=None):
    """MS1 scans with a Gaussian-ish elution of M, M+1, M+2 at the given intensity
    ratios. `extra` = (mz, profile) adds one more ion with its own elution profile."""
    step = NEUTRON_MASS / charge
    profile = [0, 5, 40, 100, 60, 15, 0]
    specs = []
    for i, h in enumerate(profile):
        peaks = [(mono_mz + k * step, h * r) for k, r in enumerate(ratios) if h * r > 0]
        if extra is not None and extra[1][i] > 0:
            peaks.append((extra[0], extra[1][i]))
        peaks.append((300.0, 10.0))
        specs.append(_ms1(2 * i + 1, 20.0 + 0.1 * i, sorted(peaks)))
    return FakeReader(specs)


class TestIsotopeTraceEvidence:
    MONO, Z = 1204.5640, 3

    def test_true_isotopes_coelute(self):
        ev = isotope_trace_evidence(_envelope_reader(self.MONO, self.Z, [1.0, 0.9, 0.5]),
                                    self.MONO, self.Z, rt_center=20.3, rt_halfwidth=1.0)
        assert ev.found
        assert ev.apex_rt == pytest.approx(20.3)
        assert ev.coelution[1] == pytest.approx(1.0)
        assert ev.coelution[2] == pytest.approx(1.0)
        assert ev.peak_intensity[-1] == 0.0
        assert not ev.m_minus_1_coelutes

    def test_coeluting_m_minus_1_is_flagged(self):
        # The "assigned" mono is really M+1 of a species one neutron lighter.
        step = NEUTRON_MASS / self.Z
        reader = _envelope_reader(self.MONO - step, self.Z, [0.6, 1.0, 0.8, 0.4])
        ev = isotope_trace_evidence(reader, self.MONO, self.Z, rt_center=20.3)
        assert ev.coelution[-1] == pytest.approx(1.0)
        assert ev.m_minus_1_coelutes

    def test_m_minus_1_that_elutes_elsewhere_is_not_flagged(self):
        step = NEUTRON_MASS / self.Z
        other = (self.MONO - step, [80, 40, 0, 0, 0, 0, 0])
        ev = isotope_trace_evidence(
            _envelope_reader(self.MONO, self.Z, [1.0, 0.9, 0.5], extra=other),
            self.MONO, self.Z, rt_center=20.3)
        assert ev.peak_intensity[-1] == 0.0
        assert ev.coelution[-1] < 0.1
        assert not ev.m_minus_1_coelutes

    def test_pattern_cosine_uses_a_supplied_distribution(self):
        from pyteomics import mass
        from mzml_utils.isotopes import isotope_distribution
        comp = mass.Composition(sequence="TPENFPSK") + mass.Composition({"Br": 1, "H": -1})
        dist = isotope_distribution(comp, n_peaks=3, charge=2)
        specs = []
        for i, h in enumerate([0, 10, 100, 30, 0]):
            peaks = [(mz, h * a) for mz, a in zip(dist.mz, dist.abundance) if h > 0]
            specs.append(_ms1(i + 1, 12.0 + 0.05 * i, peaks + [(300.0, 1.0)]))
        ev = isotope_trace_evidence(FakeReader(specs), float(dist.mz[0]), 2, rt_center=12.1,
                                    distribution=dist)
        assert ev.pattern_cosine == pytest.approx(1.0)
        # The averagine cannot express M+2 > M, so it fits the same data worse.
        avg = isotope_trace_evidence(FakeReader(specs), float(dist.mz[0]), 2, rt_center=12.1,
                                     tolerance=20.0)
        assert avg.pattern_cosine < 0.9

    def test_missing_precursor(self):
        ev = isotope_trace_evidence(_envelope_reader(800.0, 2, [1.0, 0.5]), self.MONO, self.Z,
                                    rt_center=20.3)
        assert not ev.found
        assert ev.pattern_cosine == 0.0
        assert not ev.m_minus_1_coelutes


# ---------------------------------------------------------------------------
# Peak boundaries and quality metrics (ports of Chromatograms / MsQuality /
# MetaboCoreUtils; tests/data/r_peak_reference.json holds the R values)
# ---------------------------------------------------------------------------

import json
import math
from pathlib import Path

from mzml_utils.xic import peak_boundary, peak_metrics, beta_values, PeakMetrics

_R_REF = json.loads((Path(__file__).parent / "data" / "r_peak_reference.json").read_text())


def _gauss(rt, c, s, h):
    return h * np.exp(-0.5 * ((rt - c) / s) ** 2)


class TestRReference:
    """The port reproduces the R functions run verbatim on the same traces."""

    @pytest.mark.parametrize("name", sorted(_R_REF["traces"]))
    def test_matches_r(self, name):
        tr = _R_REF["traces"][name]
        rt, y, ref = np.array(tr["rt"]), np.array(tr["intensity"]), tr["r"]
        xic = XIC(target_mz=0.0, rt=rt, intensity=y)
        left, _apex, right = peak_boundary(xic)
        assert rt[left] == ref["left_rt"] and rt[right] == ref["right_rt"]
        seg = slice(left, right + 1)
        got = dict(zip(("beta_cor_in_boundary", "beta_snr_in_boundary"), beta_values(y[seg], rt[seg])))
        got.update(zip(("beta_cor_whole", "beta_snr_whole"), beta_values(y, rt)))
        pm = peak_metrics(xic)
        got["snr"], got["prominence"] = pm.snr, pm.prominence
        assert pm.beta_cor == pytest.approx(ref["beta_cor_in_boundary"], rel=1e-12)
        for k, v in got.items():
            if ref[k] is None:  # R Inf: MAD or q10 is zero; MsQuality returns NA -> NaN here
                assert math.isnan(v), k
            else:
                assert v == pytest.approx(ref[k], rel=1e-12), k


class TestPeakBoundary:
    def test_valley_between_two_peaks(self):
        rt = np.linspace(0, 10, 101)
        y = _gauss(rt, 3, 0.4, 1e5) + _gauss(rt, 6.5, 0.4, 4e4) + 100.0
        left, apex, right = peak_boundary(XIC(0.0, rt, y))
        assert rt[apex] == pytest.approx(3.0)
        assert 4.0 < rt[right] < 6.0  # stops in the valley, not at the second peak

    def test_apex_window_selects_the_smaller_peak(self):
        rt = np.linspace(0, 10, 101)
        y = _gauss(rt, 3, 0.4, 1e5) + _gauss(rt, 6.5, 0.4, 4e4) + 100.0 + 0.01 * rt
        left, apex, right = peak_boundary(XIC(0.0, rt, y), apex_window=(6.0, 7.0))
        assert rt[apex] == pytest.approx(6.5)
        assert 4.0 < rt[left] < 6.0

    def test_contains_rt_takes_the_peak_whose_tail_holds_the_ms2(self):
        # MS2 triggered on the tail, 1 min after the apex: the containing peak is the
        # big one, not the tallest point near the MS2 time.
        rt = np.linspace(60, 72, 121)
        y = (_gauss(rt, 64, 0.8, 3e4) + _gauss(rt, 67.1, 0.25, 1.8e8)
             + np.where(rt > 67.1, 1.5e7 * np.exp(-(rt - 67.1) / 0.6), 0.0))
        xic = XIC(0.0, rt, y)
        left, apex, right = peak_boundary(xic, contains_rt=68.1, nearest_valley=True)
        assert rt[apex] == pytest.approx(67.1) and rt[left] <= 68.1 <= rt[right]
        w_left, w_apex, _ = peak_boundary(xic, apex_window=(67.6, 68.6))
        assert rt[w_apex] > 67.5            # the window picks a tail point

    def test_contains_rt_peels_to_the_smaller_peak(self):
        rt = np.linspace(0, 10, 101)
        y = _gauss(rt, 3, 0.4, 1e5) + _gauss(rt, 6.5, 0.4, 4e4) + 100.0 + 0.01 * rt
        xic = XIC(0.0, rt, y)
        _l, apex, _r = peak_boundary(xic, contains_rt=6.8)
        assert rt[apex] == pytest.approx(6.5)
        assert peak_boundary(xic, contains_rt=50.0) is None       # outside the trace
        with pytest.raises(ValueError):
            peak_boundary(xic, contains_rt=6.8, apex_window=(6, 7))

    def test_high_valley_falls_back_to_threshold_on_both_sides(self):
        rt = np.linspace(0, 10, 101)
        y = _gauss(rt, 4.6, 0.5, 1e5) + _gauss(rt, 5.6, 0.5, 9e4) + 100.0 + 0.01 * rt
        left, apex, right = peak_boundary(XIC(0.0, rt, y))
        # The valley between the two peaks is far above baseline + 10 %, so both
        # boundaries come from the 10 % threshold and enclose both peaks.
        assert rt[left] < 4.0 and rt[right] > 6.0

    def test_zero_filled_left_boundary_runs_to_trace_start_in_r_mode(self):
        rt = np.linspace(0, 10, 101)
        y = np.where(np.abs(rt - 5) < 1.2, _gauss(rt, 5, 0.35, 2e5), 0.0)
        xic = XIC(0.0, rt, y)
        left, _apex, right = peak_boundary(xic)
        assert left == 0                        # R behaviour, documented
        l2, _a2, r2 = peak_boundary(xic, nearest_valley=True)
        # both boundaries are the first zero next to the peak
        assert y[l2] == 0 and y[l2 + 1] > 0 and y[r2] == 0 and y[r2 - 1] > 0
        assert r2 == right                      # the right side is unchanged

    def test_nearest_valley_on_noisy_baseline_matches_r(self):
        tr = _R_REF["traces"]["two_peaks"]
        xic = XIC(0.0, np.array(tr["rt"]), np.array(tr["intensity"]))
        assert peak_boundary(xic) == peak_boundary(xic, nearest_valley=True)

    def test_flat_top_is_not_mistaken_for_a_valley(self):
        rt = np.arange(9, dtype=float)
        y = np.array([0, 1, 5, 9, 9, 9, 5, 1, 0], dtype=float)
        assert peak_boundary(XIC(0.0, rt, y), nearest_valley=True) == (0, 3, 8)

    def test_degenerate_inputs(self):
        assert peak_boundary(XIC(0.0, np.array([1.0, 2.0]), np.array([5.0, 3.0]))) is None
        assert peak_boundary(XIC(0.0, np.arange(5.0), np.zeros(5))) is None
        xic = XIC(0.0, np.arange(5.0), np.array([0, 1, 5, 1, 0.0]))
        assert peak_boundary(xic, apex_window=(10.0, 11.0)) is None
        with pytest.raises(ValueError):
            peak_boundary(XIC(0.0, np.arange(5.0), np.array([0, 1, np.nan, 1, 0.0])))


class TestPeakMetrics:
    def test_area_is_within_boundaries_not_whole_trace(self):
        rt = np.linspace(0, 20, 201)
        y = _gauss(rt, 5, 0.3, 1e5) + _gauss(rt, 15, 0.3, 1e5)
        xic = XIC(0.0, rt, y)
        pm = peak_metrics(xic, nearest_valley=True)
        one_peak = 1e5 * 0.3 * math.sqrt(2 * math.pi)
        assert pm.area == pytest.approx(one_peak, rel=1e-3)
        assert xic.area == pytest.approx(2 * one_peak, rel=1e-3)

    def test_asymmetry_symmetric_and_tailing(self):
        rt = np.linspace(0, 10, 1001)
        sym = peak_metrics(XIC(0.0, rt, _gauss(rt, 5, 0.3, 1e5) + 10.0))
        assert sym.asymmetry == pytest.approx(1.0, abs=0.02)
        # exponentially modified Gaussian-like tail on the right
        tail = np.where(rt >= 4, np.exp(-(rt - 4) / 0.8) * (1 - np.exp(-(rt - 4) / 0.15)), 0.0)
        pm = peak_metrics(XIC(0.0, rt, 1e5 * tail + 10.0))
        assert pm.asymmetry > 1.5

    def test_fwhm_matches_xic_property_for_the_tallest_peak(self):
        rt = np.linspace(0, 10, 101)
        xic = XIC(0.0, rt, _gauss(rt, 5, 0.4, 1e5) + 10.0)
        pm = peak_metrics(xic)
        assert pm.fwhm == pytest.approx(xic.fwhm)
        assert pm.fwhm == pytest.approx(2.3548 * 0.4, rel=0.02)
        assert pm.points_across_peak == xic.points_across_peak

    def test_beta_cor_high_for_a_peak_low_for_noise(self):
        rt = np.linspace(0, 1, 41)
        peak = np.where((rt > 0) & (rt < 1), rt ** 2 * (1 - rt) ** 4, 0.0) * 1e6 + 1.0
        assert beta_values(peak, rt)[0] > 0.99
        noise = _R_REF["traces"]["noisy_flat"]
        assert beta_values(np.array(noise["intensity"]))[0] < 0.5

    def test_no_peak_gives_empty_metrics(self):
        pm = peak_metrics(XIC(0.0, np.array([]), np.array([])))
        assert isinstance(pm, PeakMetrics) and not pm.found and math.isnan(pm.apex_rt)
        d = pm.as_dict()
        assert "width" in d and "beta_cor" in d
