#!/usr/bin/env python3
"""Self-contained unit test for mzml_utils.spectrum_cache -- NO mzML, NO network.

Builds a tiny DB from synthetic Spectrum objects (MS1 + HCD + 5-scan EThcD SA burst),
then exercises the reader: peak round-trip, precursor fields, HCD<->EThcD pairing,
preceding-MS1, find_best_ms1, extras, in-memory mode, open_spectra() dispatch, and the
copy-first _copy_to_local helper.
Run: python3 -m pytest tests/test_spectrum_cache.py   (or: python3 tests/test_spectrum_cache.py)
"""
import os
import sys
import tempfile

import numpy as np

from mzml_utils.reader import Spectrum
from mzml_utils import spectrum_cache as sc


def _mk(scan, ms_level, act, iso=0.0, mz=None, inten=None, charge=0):
    return Spectrum(
        scan_num=scan, mz=np.array(mz or [], dtype=float),
        intensity=np.array(inten or [], dtype=float), ms_level=ms_level,
        rt=float(scan), filter_string=f"FTMS + p NSI {act}", activation_type=act,
        precursor_mz=iso, precursor_charge=charge, isolation_window_target=iso,
        isolation_window_lower=1.0, isolation_window_upper=1.0,
        base_peak_mz=(mz[0] if mz else 0.0), tic=float(sum(inten or [0])),
    )


def main():
    fails = []

    def check(cond, msg):
        print(("  PASS " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    # --- peak codec round-trip ---
    mz = np.array([126.0551, 204.0867, 366.1395, 800.5], dtype=float)
    it = np.array([1000.0, 5000.0, 250.0, 99.0], dtype=float)
    mzb, itb = sc.encode_peaks(mz, it)
    dmz, dit = sc.decode_peaks(mzb, itb)
    check(np.allclose(dmz, mz, atol=1e-2) and np.allclose(dit, it, atol=1e-1),
          "peak codec round-trips (float32+zlib)")

    tmp = tempfile.mkdtemp()
    db = os.path.join(tmp, "synthetic.spectra.db")
    conn = sc.open_db(db)

    scan_meta = []
    # scan 1 MS1
    sc.insert_scan(conn, _mk(1, 1, "MS1", mz=[400.0, 500.0], inten=[10, 20]))
    scan_meta.append(dict(scan_num=1, activation_type="MS1", iso_target=0.0, precursor_mz=0.0, sa=None))
    # scan 2 HCD @800.5, precursor MS1 = scan 1
    sc.insert_scan(conn, _mk(2, 2, "HCD", iso=800.5, mz=list(mz), inten=list(it), charge=3),
                   extras=dict(precursor_spectrum_ref=1, ion_injection_time=22.5,
                               collision_energy=30.0, raw_meta={"note": "hcd"}))
    scan_meta.append(dict(scan_num=2, activation_type="HCD", iso_target=800.5, precursor_mz=800.5, sa=None))
    # scans 3-7 EThcD burst @800.5, SA ladder
    for k, (s, sa) in enumerate(zip(range(3, 8), [15, 20, 25, 30, 35])):
        sc.insert_scan(conn, _mk(s, 2, "EThcD", iso=800.5, mz=[126.05, 204.09, 500.2 + k],
                                 inten=[900, 400, 50], charge=3),
                       extras=dict(precursor_spectrum_ref=1))
        scan_meta.append(dict(scan_num=s, activation_type="EThcD", iso_target=800.5,
                              precursor_mz=800.5, sa=float(sa)))
    # scan 8 MS1, scan 9 HCD @650.2, scans 10-14 EThcD burst
    sc.insert_scan(conn, _mk(8, 1, "MS1", mz=[410.0], inten=[15]))
    scan_meta.append(dict(scan_num=8, activation_type="MS1", iso_target=0.0, precursor_mz=0.0, sa=None))
    sc.insert_scan(conn, _mk(9, 2, "HCD", iso=650.2, mz=[204.087, 650.2], inten=[800, 40], charge=2),
                   extras=dict(precursor_spectrum_ref=8))
    scan_meta.append(dict(scan_num=9, activation_type="HCD", iso_target=650.2, precursor_mz=650.2, sa=None))
    for s, sa in zip(range(10, 15), [15, 20, 25, 30, 35]):
        sc.insert_scan(conn, _mk(s, 2, "EThcD", iso=650.2, mz=[126.05, 300.1], inten=[700, 30], charge=2),
                       extras=dict(precursor_spectrum_ref=8))
        scan_meta.append(dict(scan_num=s, activation_type="EThcD", iso_target=650.2,
                              precursor_mz=650.2, sa=float(sa)))

    n_pairs = sc.finalize_pairs(conn, scan_meta)
    sc.set_meta(conn, schema_version=sc.SCHEMA_VERSION, source_mzml="/fake/synthetic.mzML",
                source_size=123, source_mtime=456, n_scans=14, n_pairs=n_pairs)
    conn.commit()
    conn.close()
    check(n_pairs == 10, f"pairing built 10 child rows (got {n_pairs})")

    # --- reader (disk, read-only) ---
    cache = sc.SpectrumCache(db, mirror_local=False)
    sp = cache.get_spectrum(2)
    check(sp is not None and sp.activation_type == "HCD" and sp.precursor_charge == 3,
          "get_spectrum returns HCD Spectrum with precursor fields")
    check(sp is not None and np.allclose(np.sort(sp.mz), np.sort(mz), atol=1e-2),
          "HCD peaks round-trip through the cache (lossless)")
    check(isinstance(sp, Spectrum), "returned object IS mzml_utils.Spectrum (drop-in)")

    kids = cache.get_children(2)
    check(kids == [3, 4, 5, 6, 7], f"get_children(HCD 2) -> 5 EThcD in order (got {kids})")
    check(cache.get_children(9) == [10, 11, 12, 13, 14], "get_children(HCD 9) -> [10..14]")
    check(cache.get_parent(5) == 2, "get_parent(EThcD 5) -> HCD 2")
    check(cache.find_preceding_ms1(2) == 1, "find_preceding_ms1(HCD 2) -> MS1 1 (via precursor ref)")
    check(cache.find_preceding_ms1(9) == 8, "find_preceding_ms1(HCD 9) -> MS1 8")

    # --- find_best_ms1 (drop-in for MzMLReader.find_best_ms1) ---
    check(cache.find_best_ms1(2, 800.5) == 1, "find_best_ms1(HCD 2) -> preceding MS1 1 (fallback)")
    check(cache.find_best_ms1(9, 650.2) == 8, "find_best_ms1(HCD 9) -> preceding MS1 8 (fallback)")

    md = cache.extract_metadata(2)
    check(md and md["activation_type"] == "HCD" and abs(md["isolation_width"] - 2.0) < 1e-6,
          "extract_metadata mirrors MzMLReader (activation + isolation_width)")
    rm = cache.raw_meta(2)
    check(rm == {"note": "hcd"}, "raw_meta catch-all round-trips")
    check(cache.meta()["n_scans"] == 14, "meta() provenance readable")
    n_ms1 = sum(1 for s in cache.iter_spectra() if s.ms_level == 1)
    check(n_ms1 == 2, f"iter_spectra sees all scans (2 MS1, got {n_ms1})")
    cache.close()

    # --- in-memory mode + open_spectra dispatch ---
    cmem = sc.SpectrumCache(db, mirror_local=False, in_memory=True)
    check(cmem.get_spectrum(3).activation_type == "EThcD", "in_memory mode serves scans")
    cmem.close()
    check(isinstance(sc.open_spectra(db, mirror_local=False), sc.SpectrumCache),
          "open_spectra('*.db') -> SpectrumCache")
    # open_spectra is also exported at package top level (the single reader entry point)
    import mzml_utils as mu
    check(isinstance(mu.open_spectra(db, mirror_local=False), sc.SpectrumCache),
          "mzml_utils.open_spectra('*.db') -> SpectrumCache")

    # --- _copy_to_local: sequential source copy preserves bytes + mtime (no mzML/network) ---
    src = os.path.join(tmp, "src.bin")
    with open(src, "wb") as fh:
        fh.write(b"spectrum-cache copy-first test payload " * 2000)
    os.utime(src, (1_000_000_000, 1_700_000_000))
    dst = os.path.join(tmp, "copied.bin")
    sc._copy_to_local(src, dst)
    check(open(dst, "rb").read() == open(src, "rb").read(), "_copy_to_local copies bytes exactly")
    check(int(os.stat(dst).st_mtime) == 1_700_000_000, "_copy_to_local preserves source mtime")

    print()
    if fails:
        print(f"RESULT: {len(fails)} FAILURE(S)")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("RESULT: ALL PASS")


def test_spectrum_cache():
    """pytest entry point."""
    main()


def _build_min_db(tmp, source_mzml, source_size=1, source_mtime=1):
    db = os.path.join(tmp, "min.spectra.db")
    conn = sc.open_db(db)
    sc.insert_scan(conn, _mk(1, 1, "MS1", mz=[400.0], inten=[10]))
    sc.set_meta(conn, schema_version=sc.SCHEMA_VERSION, source_mzml=source_mzml,
                source_size=source_size, source_mtime=source_mtime, n_scans=1, n_pairs=0)
    conn.commit()
    conn.close()
    return db


def test_ion_injection_time_round_trip(tmp_path):
    """ion_injection_time survives Spectrum -> cache -> Spectrum, and is exposed
    by extract_metadata.

    Regression guard: the field was absent from the Spectrum dataclass entirely, so
    the only way to reach it was cache-private SQL -- which breaks the contract that
    SpectrumCache and MzMLReader present an identical API. It is also capped
    per scan type on real instruments, so it must be readable alongside
    activation_type rather than separately.
    """
    db = os.path.join(str(tmp_path), "inj.spectra.db")
    conn = sc.open_db(db)

    # (a) supplied via extras, as build_cache does when parsing an mzML
    s1 = _mk(1, 1, "MS1", mz=[400.0], inten=[10])
    sc.insert_scan(conn, s1, {"ion_injection_time": 12.5})

    # (b) carried on the Spectrum with no extras -- must NOT be lost
    s2 = _mk(2, 2, "HCD", iso=500.0, mz=[204.0867], inten=[99])
    s2.ion_injection_time = 87.25
    sc.insert_scan(conn, s2)

    # (c) genuinely absent -> 0.0, never an exception
    sc.insert_scan(conn, _mk(3, 2, "EThcD", iso=600.0, mz=[366.14], inten=[5]))

    sc.set_meta(conn, schema_version=sc.SCHEMA_VERSION, source_mzml="x.mzML",
                source_size=1, source_mtime=1, n_scans=3, n_pairs=0)
    conn.commit()
    conn.close()

    cache = sc.SpectrumCache(db)
    try:
        assert cache.get_spectrum(1).ion_injection_time == 12.5
        assert cache.get_spectrum(2).ion_injection_time == 87.25
        assert cache.get_spectrum(3).ion_injection_time == 0.0
        md = cache.extract_metadata(2)
        assert md["ion_injection_time"] == 87.25
        # readable alongside the scan type it must be grouped by
        assert md["activation_type"] == "HCD"
    finally:
        cache.close()


def test_source_blind_safety():
    """P1 (safety): the cache never touches a network source by default.

    os.stat / MzMLReader against a wedged SMB share can hang uninterruptibly, so
    is_stale() and verify_cache() must refuse a /Volumes source unless explicitly allowed.
    """
    tmp = tempfile.mkdtemp()

    # network source -> is_stale refuses (None) without stat-ing; verify_cache refuses (raises)
    net_db = _build_min_db(tmp, "/Volumes/fake_share/x.mzML")
    c = sc.SpectrumCache(net_db, mirror_local=False)
    assert c.is_stale() is None, "is_stale must refuse a /Volumes source by default"
    c.close()
    try:
        sc.verify_cache(net_db)
        raise AssertionError("verify_cache must refuse a network source by default")
    except RuntimeError:
        pass

    # local source that matches -> not stale; after mutation -> stale (local path still works)
    local_src = os.path.join(tmp, "local_src.mzML")
    with open(local_src, "wb") as fh:
        fh.write(b"x" * 100)
    st = os.stat(local_src)
    loc_db = _build_min_db(tmp, local_src, source_size=st.st_size, source_mtime=int(st.st_mtime))
    c = sc.SpectrumCache(loc_db, mirror_local=False)
    assert c.is_stale() is False, "matching local source -> not stale"
    with open(local_src, "ab") as fh:
        fh.write(b"more")
    assert c.is_stale() is True, "changed local source -> stale"
    c.close()


if __name__ == "__main__":
    main()
    test_source_blind_safety()
    print("RESULT: source-blind safety PASS")
