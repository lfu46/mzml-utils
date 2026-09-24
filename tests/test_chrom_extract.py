"""Tests for chrom_extract: peak table -> chromatograms across runs."""

import math

import numpy as np
import pytest

from mzml_utils import Spectrum
from mzml_utils.chrom_extract import (extract_chromatograms, chromatogram_rows,
                                      source_key, ExtractedChromatogram)
from mzml_utils.xic import extract_xics


def _ms1(scan, rt, peaks):
    mz = np.array([p[0] for p in peaks], dtype=float)
    inten = np.array([p[1] for p in peaks], dtype=float)
    return Spectrum(scan_num=scan, mz=mz, intensity=inten, ms_level=1, rt=rt,
                    tic=float(inten.sum()))


def _ms2(scan, rt, peaks):
    mz = np.array([p[0] for p in peaks], dtype=float)
    inten = np.array([p[1] for p in peaks], dtype=float)
    return Spectrum(scan_num=scan, mz=mz, intensity=inten, ms_level=2, rt=rt,
                    precursor_mz=1000.5, activation_type="HCD")


class FakeReader:
    def __init__(self, specs):
        self._specs = specs

    def iter_spectra(self):
        yield from self._specs


def _run(shift=0.0, scale=1.0):
    """20 MS1 scans 1 min apart. Species A (1000.5000) elutes at 5, B (800.2500) at 14,
    plus a neighbour 0.02 Da from A that must stay outside a 10 ppm window."""
    specs = []
    for i in range(20):
        rt = float(i) + shift
        a = 1e5 * scale * math.exp(-0.5 * ((i - 5) / 1.2) ** 2)
        b = 5e4 * scale * math.exp(-0.5 * ((i - 14) / 1.0) ** 2)
        peaks = sorted([(1000.5000, a), (1000.5200, 7e3), (800.2500, b), (650.0, 100.0)])
        specs.append(_ms1(2 * i + 1, rt, peaks))
        specs.append(_ms2(2 * i + 2, rt + 0.5, [(204.087, 9e3)]))
    return FakeReader(specs)


TABLE = [
    {"id": "A", "mz": 1000.5, "rt_min": 2.0, "rt_max": 9.0, "run": "r1", "note": "kept"},
    {"id": "B", "mz": 800.25, "rt_min": 10.0, "rt_max": 18.0, "run": "r2", "note": "kept too"},
]


def test_source_key():
    assert source_key("/x/spectra_cache/run_01.spectra.db") == "run_01"
    assert source_key("/x/run_01.mzML") == "run_01"
    assert source_key("C:/x/run_01.mzml") == "run_01"


def test_per_row_windows_and_carry_through():
    res = extract_chromatograms(TABLE, {"r1": _run()}, metrics=True)
    assert [r.row["id"] for r in res] == ["A", "B"]
    a, b = res
    assert a.xic.rt.min() >= 2.0 and a.xic.rt.max() <= 9.0 and a.xic.n_points == 8
    assert b.xic.rt.min() >= 10.0 and b.xic.n_points == 9
    assert a.metrics.apex_rt == 5.0 and b.metrics.apex_rt == 14.0
    assert a.row == TABLE[0]                      # carried through unchanged
    assert a.xic.intensity.max() == pytest.approx(1e5)  # the 0.02 Da neighbour is excluded


def test_by_restricts_rows_to_their_run():
    srcs = {"r1": _run(), "r2": _run(shift=0.3)}
    res = extract_chromatograms(TABLE, srcs, by="run")
    assert [(r.row["id"], r.source) for r in res] == [("A", "r1"), ("B", "r2")]
    assert res[1].metrics.apex_rt == pytest.approx(14.3)


def test_without_by_every_row_from_every_run_for_overlays():
    srcs = {"r1": _run(), "r2": _run(shift=0.4, scale=0.5)}
    res = extract_chromatograms(TABLE[:1], srcs)
    assert [r.source for r in res] == ["r1", "r2"]
    assert res[1].metrics.apex_rt - res[0].metrics.apex_rt == pytest.approx(0.4)
    assert res[1].metrics.apex_intensity == pytest.approx(0.5e5)


def test_single_row_matches_extract_xics():
    reader = _run()
    res = extract_chromatograms(TABLE[:1], {"r1": reader}, tolerance=10.0, metrics=False)
    ref = extract_xics(reader, {"A": 1000.5}, tolerance=10.0, ms_level=1,
                       rt_range=(2.0, 9.0), collect_tic=False)["A"]
    np.testing.assert_array_equal(res[0].xic.rt, ref.rt)
    np.testing.assert_allclose(res[0].xic.intensity, ref.intensity, rtol=1e-12, atol=0)
    assert res[0].metrics is None


def test_unsorted_mz_is_handled():
    specs = [_ms1(1, 1.0, [(900.0, 1.0), (500.0, 7.0), (1000.5, 3.0)]),
             _ms1(2, 2.0, [(1000.5, 4.0), (100.0, 9.0)])]
    res = extract_chromatograms([{"mz": 1000.5, "rt_min": 0, "rt_max": 3}],
                                {"x": FakeReader(specs)}, metrics=False)
    np.testing.assert_allclose(res[0].xic.intensity, [3.0, 4.0])


def test_mz_min_max_override_tolerance_and_da_unit():
    row = {"mz_min": 1000.49, "mz_max": 1000.53, "rt_min": 4.0, "rt_max": 6.0}
    res = extract_chromatograms([row], {"r": _run()}, metrics=False)
    assert res[0].xic.intensity[1] == pytest.approx(1e5 + 7e3)  # both peaks inside
    assert res[0].xic.target_mz == pytest.approx(1000.51)
    res = extract_chromatograms([{"mz": 1000.5, "rt_min": 4.0, "rt_max": 6.0}], {"r": _run()},
                                tolerance=0.03, unit="Da", metrics=False)
    assert res[0].xic.intensity[1] == pytest.approx(1e5 + 7e3)


def test_apex_rt_column_measures_the_peak_containing_the_ms2():
    specs = []
    for i in range(30):
        y = 1e5 * math.exp(-0.5 * ((i - 6) / 1.0) ** 2) + 3e4 * math.exp(-0.5 * ((i - 20) / 1.0) ** 2)
        specs.append(_ms1(i + 1, float(i), [(1000.5, y + 10.0 + 0.01 * i)]))
    row = {"mz": 1000.5, "rt_min": 0.0, "rt_max": 29.0, "ms2_rt": 21.6}
    res = extract_chromatograms([row], {"r": FakeReader(specs)}, apex_rt="ms2_rt")
    assert res[0].metrics.apex_rt == 20.0
    res = extract_chromatograms([row], {"r": FakeReader(specs)})
    assert res[0].metrics.apex_rt == 6.0
    row["ms2_rt"] = float("nan")          # missing MS2 time falls back to the tallest
    res = extract_chromatograms([row], {"r": FakeReader(specs)}, apex_rt="ms2_rt")
    assert res[0].metrics.apex_rt == 6.0


def test_row_with_no_scans_gets_empty_xic():
    res = extract_chromatograms([{"mz": 1000.5, "rt_min": 100.0, "rt_max": 101.0}], {"r": _run()})
    assert res[0].xic.n_points == 0 and not res[0].metrics.found


def test_validation_errors():
    with pytest.raises(ValueError, match="rt_min"):
        extract_chromatograms([{"mz": 1.0, "rt_max": 2.0}], {"r": _run()})
    with pytest.raises(ValueError, match="NaN"):
        extract_chromatograms([{"mz": 1.0, "rt_min": float("nan"), "rt_max": 2.0}], {"r": _run()})
    with pytest.raises(ValueError, match="match no source"):
        extract_chromatograms(TABLE, {"r1": _run()}, by="run")
    with pytest.raises(ValueError, match="both"):
        extract_chromatograms([{"mz_min": 1.0, "rt_min": 0, "rt_max": 2}], {"r": _run()})
    with pytest.raises(ValueError, match="share a run name"):
        extract_chromatograms(TABLE, ["/a/x.mzML", "/b/x.mzML"])


def test_chromatogram_rows_flattens_and_warns_on_collision():
    res = extract_chromatograms(TABLE[:1], {"r1": _run()})
    rows = chromatogram_rows(res)
    assert rows[0]["id"] == "A" and rows[0]["source"] == "r1" and rows[0]["apex_rt"] == 5.0
    assert {"n_scans", "area", "beta_cor", "width"} <= set(rows[0])
    clash = [dict(TABLE[0], area=-1)]
    with pytest.warns(UserWarning, match="area"):
        rows = chromatogram_rows(extract_chromatograms(clash, {"r1": _run()}))
    assert rows[0]["area"] > 0
