"""ms1quant: the IonQuant-like top-3 isotope-cluster estimator, on synthetic runs (no files)."""
import os

import numpy as np
import pytest

from mzml_utils import ms1quant as mq
from mzml_utils import spectrum_cache as sc
from mzml_utils.reader import Spectrum

D = mq.ISOTOPE_SPACING
ANCHOR, Z = 900.0, 3


def _cluster(scale, isotopes=(0, 1, 2, 3)):
    rel = {-1: 0.0, 0: 1.0, 1: 0.9, 2: 0.5, 3: 0.2, 4: 0.07}
    return [(ANCHOR + i * D / Z, scale * rel[i]) for i in isotopes]


def _run(tmp_path, profile, filler=None, ms2_rt=10.0):
    """MS1 every 0.1 min from 9.0 to 11.0 with the cluster scaled by profile(rt); one MS2 at ms2_rt."""
    db = os.path.join(str(tmp_path), "spectra_cache", "run.spectra.db")
    os.makedirs(os.path.dirname(db), exist_ok=True)
    conn = sc.open_db(db)
    scan = 1
    for k in range(21):
        rt = round(9.0 + 0.1 * k, 3)
        peaks = [(400.0 + j, 1e4 * (j + 1)) for j in range(filler or 5)]
        s = profile(rt)
        if s:
            peaks += _cluster(s)
        peaks.sort()
        sc.insert_scan(conn, Spectrum(scan_num=scan, mz=np.array([p[0] for p in peaks]),
                                      intensity=np.array([p[1] for p in peaks]), ms_level=1, rt=rt,
                                      activation_type="MS1", filter_string="FTMS + p NSI Full ms"))
        scan += 1
        if abs(rt - ms2_rt) < 1e-9:
            sc.insert_scan(conn, Spectrum(scan_num=scan, mz=np.array([200.0]), intensity=np.array([1.0]),
                                          ms_level=2, rt=rt + 0.001, activation_type="HCD",
                                          precursor_mz=ANCHOR + D / Z, precursor_charge=Z,
                                          filter_string="FTMS + c NSI d Full ms2 900.33@hcd30.00"))
            ms2 = scan
            scan += 1
    sc.set_meta(conn, schema_version=1)
    conn.commit()
    conn.close()
    return db, ms2


def gauss(rt):
    return 1e6 * np.exp(-((rt - 10.0) / 0.12) ** 2) if abs(rt - 10.0) <= 0.3 else 0.0


def test_top3_matches_a_hand_computation(tmp_path):
    db, ms2 = _run(tmp_path, gauss)
    with sc.open_spectra(db) as r:
        q = mq.quantify_precursors(r, [dict(ms2_scan=ms2, charge=Z, anchor_mz=ANCHOR)])[0]
    per_scan = sorted(gauss(round(9.0 + 0.1 * k, 3)) * (1.0 + 0.9 + 0.5 + 0.2) for k in range(21))
    assert q["top"] == pytest.approx(sum(per_scan[-3:]), rel=1e-6)
    assert q["apex"] == pytest.approx(per_scan[-1], rel=1e-6)
    assert q["n_scans"] == sum(v > 0 for v in per_scan) == 5 and q["apex_rt"] == pytest.approx(10.0)
    assert q["cap_censored"] is False and q["apex_floor_ratio"] == float("inf")


def test_anchor_defaults_to_the_selected_ion(tmp_path):
    db, ms2 = _run(tmp_path, gauss)
    with sc.open_spectra(db) as r:
        q = mq.quantify_precursors(r, [dict(ms2_scan=ms2, charge=Z)])[0]
    assert q["anchor_source"] == "selected_ion" and q["anchor_mz"] == pytest.approx(ANCHOR + D / Z)
    assert q["top"] > 0                     # iso_lo = -1 still reaches the monoisotope


def test_min_scans_and_min_isotopes(tmp_path):
    narrow = lambda rt: 1e6 if abs(rt - 10.0) < 0.15 else 0.0     # 3 scans: 9.9, 10.0, 10.1
    db, ms2 = _run(tmp_path, narrow)
    with sc.open_spectra(db) as r:
        idx = mq.Ms1Index(r)
    rtc = idx.scan_rt[ms2]
    assert mq.precursor_ms1(idx, ANCHOR, Z, rtc).n_scans == 3
    assert mq.precursor_ms1(idx, ANCHOR, Z, rtc, mq.Ms1QuantParams(min_scans=4)).top == 0.0
    # one isotope only is not a cluster
    assert mq.precursor_ms1(idx, ANCHOR, Z, rtc, mq.Ms1QuantParams(min_isotopes=5)).top == 0.0
    # a wrong charge finds nothing (z = 5: no spacing D/5 multiple below +4 lands on a z = 3 isotope)
    assert mq.precursor_ms1(idx, ANCHOR, 5, rtc).top == 0.0


def test_nearest_peak_rule():
    mz = np.array([899.9950, 900.0000, 900.0060])
    assert mq._nearest_peak(mz, 900.0010, 900 * 10e-6) == 1
    assert mq._nearest_peak(mz, 900.0030, 900 * 10e-6) == 2          # tie -> the higher m/z
    assert mq._nearest_peak(mz, 900.0200, 900 * 10e-6) == -1


def test_censoring_on_a_capped_view(tmp_path):
    weak = lambda rt: 3e4 if abs(rt - 10.0) <= 0.3 else 0.0           # near the filler peaks
    db, ms2 = _run(tmp_path, weak, filler=40)
    with sc.open_spectra(db, ms1_cap=30) as r:
        assert sc.ms1_cap_status(r)["capped"]
        q = mq.quantify_precursors(r, [dict(ms2_scan=ms2, charge=Z, anchor_mz=ANCHOR)])[0]
    assert q["cap_censored"] is True


def test_params_from_workflow(tmp_path):
    p = tmp_path / "fragpipe.workflow"
    p.write_text("# comment\nionquant.mztol=12\nionquant.rttol=0.5\nionquant.minscans=2\n"
                 "ionquant.minisotopes=1\nmsfragger.precursor_true_tolerance=20\n")
    prm = mq.Ms1QuantParams.from_workflow(str(p))
    assert (prm.mztol_ppm, prm.rttol_min, prm.min_scans, prm.min_isotopes) == (12.0, 0.5, 2, 1)
    assert prm.topn == 3
