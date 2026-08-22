"""Tests for the xic module."""

import numpy as np
import pytest

from mzml_utils import Spectrum
from mzml_utils.xic import extract_xic, extract_xics, XIC, ChromatogramSet


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
