"""patch_ms1 / ms1_cap_status / the MS1 export format -- synthetic, no RAW needed.

A last test runs against a real Thermo RAW when MZML_UTILS_TEST_RAW names one whose spectrum
cache exists beside its mzML (it never writes: dry run only), and skips otherwise.
"""
import hashlib
import os
import shutil
import sqlite3
import tempfile

import numpy as np
import pytest

from mzml_utils import spectrum_cache as sc
from mzml_utils import thermo_raw as tr
from mzml_utils.reader import Spectrum

CAP = 60


def _raw_scan(scan, rt, n, seed):
    rng = np.random.default_rng(seed)
    mz = np.sort(rng.uniform(400, 1800, n))
    it = rng.uniform(1e3, 1e7, n).astype(np.float32).astype(np.float64)   # float32-exact
    return tr.RawMs1Scan(scan, rt, 1, mz, it)


def _capped(x, cap=CAP):
    keep = np.sort(np.argsort(x.intensity)[::-1][:cap])
    return x.mz[keep], x.intensity[keep]


def _spec(scan, level, mz, it, rt):
    return Spectrum(scan_num=scan, mz=np.asarray(mz, float), intensity=np.asarray(it, float),
                    ms_level=level, rt=rt, filter_string="FTMS + p NSI Full ms" if level == 1 else "FTMS + c NSI d Full ms2 800.0@hcd30.00",
                    activation_type="MS1" if level == 1 else "HCD", precursor_mz=0.0 if level == 1 else 800.0,
                    precursor_charge=0 if level == 1 else 3)


def _make_cache(tmp, raw, n_ms2=3, tamper=None):
    """A cache whose MS1 rows hold the RAW scans capped at CAP; MS2 rows in between."""
    db = os.path.join(tmp, "spectra_cache", "run.spectra.db")
    os.makedirs(os.path.dirname(db), exist_ok=True)
    conn = sc.open_db(db)
    scan = 1
    for x in raw.values():
        mz, it = _capped(x)
        if tamper:
            mz, it = tamper(x.scan, mz, it)
        sc.insert_scan(conn, _spec(x.scan, 1, mz, it, x.rt))
        for k in range(1, n_ms2 + 1):
            sc.insert_scan(conn, _spec(x.scan + k, 2, [200.0 + k, 300.0], [10.0, 20.0], x.rt + 0.001 * k))
    sc.set_meta(conn, schema_version=1)
    conn.commit()
    conn.close()
    return db


def _md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


@pytest.fixture
def raw():
    return {s: _raw_scan(s, s / 100.0, n, s) for s, n in ((1, 150), (5, 40), (9, 300), (13, CAP))}


def test_ms1x_round_trip(tmp_path, raw):
    p = tmp_path / "x.ms1x"
    tr.write_ms1x(str(p), raw.values())
    back = tr.read_ms1x(str(p))
    assert sorted(back) == sorted(raw)
    for s in raw:
        assert np.array_equal(back[s].mz, raw[s].mz) and np.array_equal(back[s].intensity, raw[s].intensity)
        assert back[s].rt == raw[s].rt and back[s].source == 1


def test_ms1x_truncated_raises(tmp_path, raw):
    p = tmp_path / "x.ms1x"
    tr.write_ms1x(str(p), raw.values())
    data = open(p, "rb").read()
    open(p, "wb").write(data[:-20])
    with pytest.raises(ValueError):
        tr.read_ms1x(str(p))


def test_check_accepts_a_most_intense_cap(raw):
    x = raw[1]
    mz, it = _capped(x)
    assert sc.check_ms1_against_raw(mz.astype(np.float32), it.astype(np.float32), x.mz, x.intensity) is None
    # an uncapped scan (RAW <= cap) is identical
    y = raw[5]
    assert sc.check_ms1_against_raw(y.mz, y.intensity, y.mz, y.intensity) is None


def test_check_rejects_other_differences(raw):
    x = raw[1]
    mz, it = _capped(x)
    shifted = mz.copy(); shifted[3] += 0.01
    assert "no RAW peak" in sc.check_ms1_against_raw(shifted, it, x.mz, x.intensity)
    scaled = it.copy(); scaled[0] *= 1.01
    assert "intensity differs" in sc.check_ms1_against_raw(mz, scaled, x.mz, x.intensity)
    # keep the weakest peaks instead of the strongest: not a most-intense cap
    weak = np.sort(np.argsort(x.intensity)[:CAP])
    assert "not a most-intense cap" in sc.check_ms1_against_raw(x.mz[weak], x.intensity[weak], x.mz, x.intensity)
    assert "empty" in sc.check_ms1_against_raw([], [], x.mz, x.intensity)


def test_patch_replaces_ms1_and_leaves_ms2(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw)
    with sc.SpectrumCache(db, mirror_local=False) as c:
        before_ms2 = {s: c.get_spectrum(s) for s in (2, 3, 4)}
        assert sc.ms1_cap_status(c)["capped"] is True and sc.ms1_cap_status(c)["cap"] == CAP
    rep = sc.patch_ms1(db, export=raw)
    assert rep["status"] == "patched"
    assert rep["ms1_truncated_before"] == 2 and rep["cap_undone"] == CAP
    assert rep["ms1_peaks_after"] == sum(len(x.mz) for x in raw.values())
    with sc.SpectrumCache(db, mirror_local=False) as c:
        for s, x in raw.items():
            sp = c.get_spectrum(s)
            assert sp.n_peaks == len(x.mz)
            assert np.allclose(sp.mz, x.mz, rtol=1e-7) and np.allclose(sp.intensity, x.intensity, rtol=1e-6)
        for s, old in before_ms2.items():
            new = c.get_spectrum(s)
            assert np.array_equal(new.mz, old.mz) and np.array_equal(new.intensity, old.intensity)
        st = sc.ms1_cap_status(c)
        assert st["patched"] is True and st["capped"] is False
        assert c.meta()["ms1_patch"]["ms2_digest"] == rep["ms2_digest"]


def test_patch_is_refused_twice_without_force(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw)
    sc.patch_ms1(db, export=raw)
    with pytest.raises(sc.PatchError):
        sc.patch_ms1(db, export=raw)
    assert sc.patch_ms1(db, export=raw, force=True)["status"] == "patched"


@pytest.mark.parametrize("case", ["rt", "scan_set", "foreign_peak", "dry_run"])
def test_failed_or_dry_patch_writes_nothing(tmp_path, raw, case):
    tamper = None
    exp = dict(raw)
    if case == "rt":
        exp[9] = tr.RawMs1Scan(9, raw[9].rt + 0.01, 1, raw[9].mz, raw[9].intensity)
    elif case == "scan_set":
        exp.pop(13)
    elif case == "foreign_peak":
        def tamper(scan, mz, it):
            if scan == 9:
                mz = mz.copy(); mz[0] = 401.23456
                order = np.argsort(mz); mz, it = mz[order], it[order]
            return mz, it
    db = _make_cache(str(tmp_path), raw, tamper=tamper)
    h = _md5(db)
    if case == "dry_run":
        rep = sc.patch_ms1(db, export=exp, dry_run=True)
        assert rep["status"].startswith("checked") and rep["cap_undone"] == CAP
    else:
        with pytest.raises(sc.PatchError):
            sc.patch_ms1(db, export=exp)
    assert _md5(db) == h
    assert not [f for f in os.listdir(os.path.dirname(db)) if f.endswith(".tmp")]


def test_profile_scan_is_refused(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw)
    exp = dict(raw)
    exp[5] = tr.RawMs1Scan(5, raw[5].rt, 3, raw[5].mz, raw[5].intensity)
    with pytest.raises(sc.PatchError, match="no centroids"):
        sc.patch_ms1(db, export=exp)


def test_default_raw_is_found_beside_the_mzml(tmp_path):
    (tmp_path / "run_1.raw").write_bytes(b"x")
    assert sc._default_raw_for(str(tmp_path / "run_1_calibrated.mzML")) == str(tmp_path / "run_1.raw")
    db = tmp_path / "spectra_cache" / "run_1.spectra.db"
    assert sc._default_raw_for(str(db)) == str(tmp_path / "run_1.raw")


def test_mzml_header_cap_is_read(tmp_path):
    p = tmp_path / "x.mzML"
    p.write_text('<mzML><dataProcessingList><dataProcessing id="d"><processingMethod order="2">'
                 '<userParam name="most intense count (excluding ties at the threshold)" value="1000"/>'
                 '</processingMethod></dataProcessing></dataProcessingList></mzML>')

    class FakeReader:
        path = str(p)
    st = sc.ms1_cap_status(FakeReader())
    assert st["capped"] is True and st["cap"] == 1000 and st["source"] == "mzml-header"


RAW = os.environ.get("MZML_UTILS_TEST_RAW")


@pytest.mark.skipif(not RAW or not tr.raw_reader_available(), reason="set MZML_UTILS_TEST_RAW to a Thermo RAW with a cache beside its mzML")
def test_real_raw_export_matches_its_capped_cache():
    db = sc.cache_path_for(os.path.splitext(RAW)[0] + ".mzML")
    assert os.path.isfile(db), db
    rep = sc.patch_ms1(db, RAW, dry_run=True)
    assert rep["n_ms1"] > 0 and rep["ms1_peaks_after"] >= rep["ms1_peaks_before"]


def _shift_ppm(mz):
    return -4.0 + 4.5 * np.tanh((mz - 1100.0) / 150.0)        # smooth, like MSFragger's curve


def _calibrated_tamper(rng_seed=None):
    def tamper(scan, mz, it):
        if rng_seed is None:
            return mz * (1 + _shift_ppm(mz) / 1e6), it
        rng = np.random.default_rng(rng_seed + scan)
        return mz * (1 + rng.uniform(-8, 8, len(mz)) / 1e6), it    # not a function of m/z
    return tamper


def _set_source(db, name):
    conn = sqlite3.connect(db)
    sc.set_meta(conn, source_mzml=f"/data/{name}")
    conn.commit(); conn.close()


def test_calibrated_cache_is_uncapped_on_its_own_scale(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw, tamper=_calibrated_tamper())
    _set_source(db, "run_calibrated.mzML")
    with sc.SpectrumCache(db, mirror_local=False) as c:
        kept = {s: c.get_spectrum(s) for s in raw}
    with pytest.raises(sc.PatchError, match="calibrated"):
        sc.patch_ms1(db, export=raw, mz_mode="exact")
    rep = sc.patch_ms1(db, export=raw)                  # auto -> calibrated by the source name
    assert rep["mz_mode"] == "calibrated" and rep["heldout_p99_ppm_worst_scan"] < 0.5
    with sc.SpectrumCache(db, mirror_local=False) as c:
        for s, x in raw.items():
            sp = c.get_spectrum(s)
            assert sp.n_peaks == len(x.mz)
            truth = x.mz * (1 + _shift_ppm(x.mz) / 1e6)
            err = np.abs(np.sort(sp.mz.astype(float)) - np.sort(truth)) / truth * 1e6
            assert err.max() < 0.5                       # float32 storage + interpolation
            assert set(kept[s].mz.tolist()) <= set(sp.mz.tolist())   # kept peaks unchanged


def test_calibrated_mode_refuses_a_shift_that_is_not_smooth(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw, tamper=_calibrated_tamper(rng_seed=7))
    h = _md5(db)
    with pytest.raises(sc.PatchError, match="not smooth|reorders"):
        sc.patch_ms1(db, export=raw, mz_mode="calibrated")
    assert _md5(db) == h


def test_cap_peaks_excludes_ties_at_the_threshold():
    mz = np.arange(6, dtype=float)
    it = np.array([5.0, 9.0, 3.0, 3.0, 7.0, 1.0])
    assert cap_list(mz, it, 3) == [1.0, 4.0, 0.0]          # 9, 7, 5
    assert cap_list(mz, it, 4) == [1.0, 4.0, 0.0]          # the two 3s tie at the 4th place: both out
    assert cap_list(mz, it, 5) == [1.0, 4.0, 0.0, 2.0, 3.0]
    assert len(cap_list(mz, it, 10)) == 6                  # fewer peaks than the cap: all kept


def cap_list(mz, it, n):
    m, i = sc.cap_peaks(mz, it, n)
    return [float(x) for x in m[np.argsort(-i, kind="stable")]]


def test_recapped_view_reproduces_the_capped_cache(tmp_path, raw):
    db = _make_cache(str(tmp_path), raw)
    with sc.SpectrumCache(db, mirror_local=False) as c:
        before = {s: c.get_spectrum(s) for s in raw}
    rep = sc.patch_ms1(db, export=raw)
    assert rep["recap_reproduces_cache"] is True
    with sc.SpectrumCache(db, mirror_local=False, ms1_cap=CAP) as c:
        for s, old in before.items():
            new = c.get_spectrum(s)
            assert np.array_equal(new.mz, old.mz) and np.array_equal(new.intensity, old.intensity)
        assert sc.ms1_cap_status(c)["capped"] is True
    with sc.open_spectra(db, ms1_cap=CAP) as c:
        assert all(sp.n_peaks <= CAP for sp in c.iter_spectra(ms_level=1))


def test_find_pwiz_dir_copies_the_dlls_locally(tmp_path, monkeypatch):
    src = tmp_path / "bundle"
    src.mkdir()
    for name in tr.DLLS:
        (src / name).write_bytes(b"dll")
    monkeypatch.setenv("MZML_UTILS_CACHE", str(tmp_path / "cache"))
    monkeypatch.setenv("THERMO_PWIZ_DIR", os.pathsep.join([str(tmp_path / "missing"), str(src)]))
    monkeypatch.setattr(tr, "LOCAL_DLL_DIRS", ())
    found = tr.find_pwiz_dir()
    assert found == tmp_path / "cache" / "thermo" / "dlls"
    assert all((found / n).read_bytes() == b"dll" for n in tr.DLLS)
    monkeypatch.delenv("THERMO_PWIZ_DIR")
    assert tr.find_pwiz_dir() == found                  # the local copy now suffices
    shutil.rmtree(found)
    assert tr.find_pwiz_dir(required=False) is None
    with pytest.raises(FileNotFoundError):
        tr.find_pwiz_dir()


def test_uncut_scan_is_not_held_to_the_smoothness_check():
    rng = np.random.default_rng(3)
    rmz = np.sort(rng.uniform(400, 1800, 40)); rit = rng.uniform(1e3, 1e6, 40).astype(np.float32).astype(float)
    cmz = rmz * (1 + rng.uniform(-8, 8, 40) / 1e6)      # rough shift, but every peak kept
    mz, it, info = sc.merge_calibrated_ms1(cmz.astype(np.float32), rit.astype(np.float32), rmz, rit)
    assert info["n_added"] == 0 and len(mz) == 40
