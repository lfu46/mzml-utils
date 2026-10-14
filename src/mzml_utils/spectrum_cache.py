#!/usr/bin/env python3
"""
spectrum-cache: build-once, read-fast local cache of everything spectrum annotation needs.

Reads an mzML ONE time (the only slow / network step) and writes a compact, float32
(near-lossless), indexed SQLite cache. All later annotation reads go through ``SpectrumCache``, which
returns ``mzml_utils.Spectrum`` objects -- a drop-in for ``mzml_utils.MzMLReader`` -- so
both the ``nglyco-analysis`` and ``oglyco-analysis`` skills work unchanged via
``open_spectra()``.

Design (see PLAN):
  * RAW peaks stored as float32 + zlib (near-lossless: not bit-exact if the source was
    float64, but <~0.2 ppm m/z error -- negligible for glyco matching). N-glyco matches
    raw (un-deisotoped) spectra and O-glyco deisotopes on the fly -- a filtered/deisotoped
    cache would break N-glyco.
  * Captures MORE than the Spectrum dataclass: FAIMS CV,
    collision energy, scan window, precursor->master-MS1 ref, polarity, centroid flag,
    plus a zlib(JSON) catch-all of every remaining scalar param (nothing lost).
  * Precomputes HCD<->EThcD (xN SA burst) pairing (O-glyco 3-panel; N-glyco ignores it).
  * One DB per mzML. Build LOCAL then copy to network (SQLite writes over SMB are unsafe).
    SpectrumCache mirrors a network DB to local disk (or :memory:) on first open so reads
    never pay SMB random-access latency.

Glyco-type-agnostic: it stores spectra + pairing only; all N/O glyco logic stays in the
consumer skills.
"""
from __future__ import annotations

import argparse
import bisect
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import zlib
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

# Reuse mzml_utils' own extraction + data model (never reimplement, never modify it).
from .reader import Spectrum, classify_activation, _extract_spectrum

SCHEMA_VERSION = 1
_HCD_RE = re.compile(r"@hcd([\d.]+)")
_ETD_RE = re.compile(r"@etd([\d.]+)")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS scans (
    scan_num              INTEGER PRIMARY KEY,
    ms_level              INTEGER,
    activation_type       TEXT,
    rt                    REAL,
    filter_string         TEXT,
    precursor_mz          REAL,
    precursor_charge      INTEGER,
    precursor_intensity   REAL,
    precursor_spectrum_ref INTEGER,   -- master MS1 scan the precursor came from
    iso_target            REAL,
    iso_lower             REAL,
    iso_upper             REAL,
    ion_injection_time    REAL,
    collision_energy      REAL,
    compensation_voltage  REAL,
    scan_window_lo        REAL,
    scan_window_hi        REAL,
    polarity              TEXT,
    centroid              INTEGER,
    tic                   REAL,
    base_peak_mz          REAL,
    base_peak_intensity   REAL,
    n_peaks               INTEGER,
    mz_blob               BLOB,        -- zlib(float32[])  RAW peaks (float32; not bit-exact)
    int_blob              BLOB,        -- zlib(float32[])  RAW peaks (float32; not bit-exact)
    raw_meta_z            BLOB         -- zlib(JSON) catch-all of remaining scalar params
);
CREATE INDEX IF NOT EXISTS idx_scans_mslevel ON scans(ms_level);
CREATE TABLE IF NOT EXISTS pairs (
    hcd_scan         INTEGER,
    child_scan       INTEGER,
    child_activation TEXT,
    sa_energy        REAL,
    order_in_burst   INTEGER
);
CREATE INDEX IF NOT EXISTS idx_pairs_hcd   ON pairs(hcd_scan);
CREATE INDEX IF NOT EXISTS idx_pairs_child ON pairs(child_scan);
"""

_SCAN_COLS = [
    "scan_num", "ms_level", "activation_type", "rt", "filter_string",
    "precursor_mz", "precursor_charge", "precursor_intensity", "precursor_spectrum_ref",
    "iso_target", "iso_lower", "iso_upper", "ion_injection_time", "collision_energy",
    "compensation_voltage", "scan_window_lo", "scan_window_hi", "polarity", "centroid",
    "tic", "base_peak_mz", "base_peak_intensity", "n_peaks", "mz_blob", "int_blob",
    "raw_meta_z",
]


# --------------------------------------------------------------------------- #
# peak / metadata (de)serialization -- small, testable, no I/O
# --------------------------------------------------------------------------- #
def encode_peaks(mz, intensity) -> Tuple[bytes, bytes]:
    mz = np.asarray(mz, dtype=np.float32)
    it = np.asarray(intensity, dtype=np.float32)
    return zlib.compress(mz.tobytes(), 6), zlib.compress(it.tobytes(), 6)


def decode_peaks(mz_blob, int_blob) -> Tuple[np.ndarray, np.ndarray]:
    mz = np.frombuffer(zlib.decompress(mz_blob), dtype=np.float32) if mz_blob else np.array([], dtype=np.float32)
    it = np.frombuffer(zlib.decompress(int_blob), dtype=np.float32) if int_blob else np.array([], dtype=np.float32)
    return mz, it


def _scan_num_from_ref(ref: Optional[str]) -> Optional[int]:
    if ref and "scan=" in ref:
        try:
            return int(ref.split("scan=")[-1].split()[0])
        except ValueError:
            return None
    return None


def _jsonable(obj, _depth=0):
    """Recursively keep only JSON-serializable scalars; drop big arrays (peaks)."""
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist() if obj.size <= 16 else None
    if _depth > 6:
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): _jsonable(v, _depth + 1) for k, v in obj.items()
                if k not in ("m/z array", "intensity array")}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v, _depth + 1) for v in obj]
    return str(obj)


def extract_extras(spec: dict) -> dict:
    """Pull the mzML fields the Spectrum dataclass does NOT carry.

    ``ion_injection_time`` is kept here too: it is stored in its own column and is
    also surfaced on Spectrum, so the cache and MzMLReader agree.
    """
    e: Dict[str, object] = {}
    scan = (spec.get("scanList", {}).get("scan") or [{}])[0]
    e["ion_injection_time"] = scan.get("ion injection time")
    e["compensation_voltage"] = scan.get("FAIMS compensation voltage",
                                         scan.get("compensation voltage"))
    sw = (scan.get("scanWindowList", {}).get("scanWindow") or [{}])[0]
    e["scan_window_lo"] = sw.get("scan window lower limit")
    e["scan_window_hi"] = sw.get("scan window upper limit")
    prec_list = spec.get("precursorList", {}).get("precursor") or []
    if prec_list:
        prec = prec_list[0]
        e["precursor_spectrum_ref"] = _scan_num_from_ref(prec.get("spectrumRef"))
        e["collision_energy"] = (prec.get("activation", {}) or {}).get("collision energy")
    e["polarity"] = ("positive" if "positive scan" in spec
                     else "negative" if "negative scan" in spec else None)
    e["centroid"] = 1 if "centroid spectrum" in spec else 0
    e["raw_meta"] = _jsonable(spec)
    return e


def _row_from(spectrum: Spectrum, extras: Optional[dict] = None) -> tuple:
    extras = extras or {}
    mz_blob, int_blob = encode_peaks(spectrum.mz, spectrum.intensity)
    raw_meta = extras.get("raw_meta")
    raw_z = zlib.compress(json.dumps(raw_meta).encode(), 6) if raw_meta is not None else None

    def num(v):
        try:
            return float(v) if v is not None else None
        except (TypeError, ValueError):
            return None

    # Spectrum now carries ion_injection_time too. Prefer the freshly parsed `extras`
    # value, but fall back to the dataclass so a Spectrum inserted without extras does
    # not silently lose it.
    inj = extras.get("ion_injection_time")
    if inj is None:
        inj = spectrum.ion_injection_time or None

    return (
        int(spectrum.scan_num), int(spectrum.ms_level), spectrum.activation_type,
        num(spectrum.rt), spectrum.filter_string,
        num(spectrum.precursor_mz), int(spectrum.precursor_charge or 0),
        num(spectrum.precursor_intensity), extras.get("precursor_spectrum_ref"),
        num(spectrum.isolation_window_target), num(spectrum.isolation_window_lower),
        num(spectrum.isolation_window_upper), num(inj),
        num(extras.get("collision_energy")), num(extras.get("compensation_voltage")),
        num(extras.get("scan_window_lo")), num(extras.get("scan_window_hi")),
        extras.get("polarity"), int(extras.get("centroid") or 0),
        num(spectrum.tic), num(spectrum.base_peak_mz), num(spectrum.base_peak_intensity),
        int(spectrum.n_peaks), mz_blob, int_blob, raw_z,
    )


def _spectrum_from_row(row: sqlite3.Row) -> Spectrum:
    mz, it = decode_peaks(row["mz_blob"], row["int_blob"])
    return Spectrum(
        scan_num=row["scan_num"], mz=mz, intensity=it, ms_level=row["ms_level"],
        rt=row["rt"] or 0.0, filter_string=row["filter_string"] or "",
        precursor_mz=row["precursor_mz"] or 0.0, precursor_charge=row["precursor_charge"] or 0,
        precursor_intensity=row["precursor_intensity"] or 0.0,
        activation_type=row["activation_type"] or "", tic=row["tic"] or 0.0,
        base_peak_mz=row["base_peak_mz"] or 0.0, base_peak_intensity=row["base_peak_intensity"] or 0.0,
        isolation_window_target=row["iso_target"] or 0.0,
        isolation_window_lower=row["iso_lower"] or 0.0,
        isolation_window_upper=row["iso_upper"] or 0.0,
        ion_injection_time=row["ion_injection_time"] or 0.0,
    )


# --------------------------------------------------------------------------- #
# low-level DB helpers (shared by build_cache and the unit test)
# --------------------------------------------------------------------------- #
def open_db(path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.executescript(_SCHEMA)
    return conn


def insert_scan(conn: sqlite3.Connection, spectrum: Spectrum, extras: Optional[dict] = None) -> None:
    ph = ",".join("?" * len(_SCAN_COLS))
    conn.execute(f"INSERT OR REPLACE INTO scans ({','.join(_SCAN_COLS)}) VALUES ({ph})",
                 _row_from(spectrum, extras))


def _iso_of(scan: dict) -> float:
    v = scan.get("iso_target") or 0.0
    return v if v else (scan.get("precursor_mz") or 0.0)


def compute_pairs(scan_meta: List[dict]) -> List[tuple]:
    """HCD<->EThcD burst pairing (ported from the validated verify_v2 logic).

    scan_meta: list of dicts with scan_num, activation_type, iso_target, precursor_mz,
    sa (float | None), sorted by scan_num. Groups EThcD/ETD into bursts by isolation m/z
    + SA-ladder reset (pd children are QUEUED, not scan-adjacent), then matches each burst
    to the nearest preceding same-isolation HCD.
    """
    scan_meta = sorted(scan_meta, key=lambda s: s["scan_num"])
    children = [s for s in scan_meta if s["activation_type"] in ("EThcD", "ETD")]
    bursts: List[dict] = []
    cur = None
    for r in children:
        iso = _iso_of(r)
        sa = r.get("sa")
        start_new = (cur is None or abs(iso - cur["iso"]) > 0.002
                     or (cur["sa"] and sa is not None and sa <= cur["sa"][-1]))
        if start_new:
            cur = dict(iso=iso, first=r["scan_num"], scans=[r["scan_num"]], sa=[sa])
            bursts.append(cur)
        else:
            cur["scans"].append(r["scan_num"])
            cur["sa"].append(sa)

    hcd = [s for s in scan_meta if s["activation_type"] == "HCD"]
    hcd_scans = [h["scan_num"] for h in hcd]
    pairs: List[tuple] = []
    for b in bursts:
        i = bisect.bisect_left(hcd_scans, b["first"]) - 1
        found = None
        steps = 0
        while i >= 0 and steps < 80:
            if abs(_iso_of(hcd[i]) - b["iso"]) <= 0.002:
                found = hcd[i]
                break
            i -= 1
            steps += 1
        if found:
            for order, (cs, sa) in enumerate(zip(b["scans"], b["sa"])):
                pairs.append((found["scan_num"], cs, "EThcD", sa, order))
    return pairs


def finalize_pairs(conn: sqlite3.Connection, scan_meta: List[dict]) -> int:
    pairs = compute_pairs(scan_meta)
    conn.execute("DELETE FROM pairs")
    conn.executemany(
        "INSERT INTO pairs (hcd_scan, child_scan, child_activation, sa_energy, order_in_burst) "
        "VALUES (?,?,?,?,?)", pairs)
    return len(pairs)


def set_meta(conn: sqlite3.Connection, **kw) -> None:
    conn.executemany("INSERT OR REPLACE INTO meta (key, value) VALUES (?,?)",
                     [(k, json.dumps(v)) for k, v in kw.items()])


def _is_network(path: str) -> bool:
    return os.path.abspath(path).startswith("/Volumes/")


def cache_path_for(mzml_path: str, subdir: str = "spectra_cache") -> str:
    """Deterministic co-located cache location: <mzml_dir>/spectra_cache/<stem>.spectra.db.

    Keeps caches grouped in one subfolder next to the mzMLs (not loose among the raws).
    """
    d = os.path.dirname(os.path.abspath(mzml_path))
    stem = os.path.splitext(os.path.basename(mzml_path))[0]
    return os.path.join(d, subdir, stem + ".spectra.db")


# --------------------------------------------------------------------------- #
# build (the ONE mzML read)
# --------------------------------------------------------------------------- #
def _build_tmpdir() -> str:
    d = os.path.join(tempfile.gettempdir(), "spectrum_cache_build")
    os.makedirs(d, exist_ok=True)
    return d


def _copy_to_local(src: str, dst: str, log=None) -> None:
    """Sequential copy src -> dst, preserving mtime. Prefer rsync (SMB-friendly and
    resumable via --partial); fall back to shutil.copy2. A plain sequential read is what
    a flaky SMB share tolerates -- unlike pyteomics' index-seek, which stalls in
    uninterruptible I/O. Raises on a size mismatch (a still-unreadable source region)."""
    t0 = time.time()
    used = "shutil"
    rsync = shutil.which("rsync")
    if rsync:
        rc = subprocess.run([rsync, "-t", "--partial", "--inplace", src, dst]).returncode
        used = "rsync" if rc == 0 else "shutil"
        if rc != 0 and log:
            log(f"rsync rc={rc}; falling back to shutil.copy2")
    if used == "shutil":
        shutil.copy2(src, dst)
    ssz, dsz = os.path.getsize(src), os.path.getsize(dst)
    if ssz != dsz:
        raise IOError(f"local copy size mismatch (src {ssz} vs dst {dsz}) -- source region unreadable?")
    if log:
        dt = time.time() - t0
        log(f"copied source to local ({used}): {dsz} B in {dt:.0f}s ({dsz/1e6/max(1e-6, dt):.1f} MB/s)")


def build_cache(mzml_path: str, db_path: str, commit_every: int = 2000,
                progress=None, keep_local_temp: bool = False,
                local_copy: Optional[bool] = None, log=None) -> dict:
    """Read an mzML once and write the cache.

    COPY-FIRST (default for a network source): a source mzML on a share (/Volumes/...)
    is first copied to local disk with a plain sequential read (rsync), and the cache is
    built from that local copy. pyteomics' ``use_index`` seeks to the file tail to read
    the mzML index -- over a flaky SMB share that pattern stalls in uninterruptible I/O
    (and hard-errors if the tail is corrupt), whereas a sequential copy is SMB-friendly,
    resumable, and doubles as a full-file integrity check. ``local_copy`` overrides:
    None = auto (True iff the source is on /Volumes), True = always copy first,
    False = read directly over the network (old behaviour).

    If db_path is on the network the DB is BUILT on local disk (SQLite writes over SMB are
    unsafe) then copied to the network. ALL local temp files -- the copied mzML and the
    local DB build file -- are removed on exit, on success OR failure, unless
    keep_local_temp=True. The final DB records the ORIGINAL (network) mzML as its source,
    so staleness detection stays correct regardless of the local copy.
    """
    from pyteomics import mzml  # only used here, in the one-time builder

    def _log(m):
        if log:
            log(m)

    orig_mzml = os.path.abspath(mzml_path)
    if local_copy is None:
        local_copy = _is_network(orig_mzml)
    dst_on_net = _is_network(db_path)
    os.makedirs(os.path.dirname(os.path.abspath(db_path)) or ".", exist_ok=True)

    tmpdir = _build_tmpdir()
    work = os.path.join(tmpdir, os.path.basename(db_path)) if dst_on_net else db_path
    tmp = work + ".tmp"
    local_mzml = None          # local copy of the source, to clean up
    conn = None
    t0 = time.time()
    try:
        # ---- source side: optional copy-to-local (the SMB-safe read) ----
        read_path = orig_mzml
        if local_copy:
            local_mzml = os.path.join(tmpdir, os.path.splitext(os.path.basename(orig_mzml))[0] + ".mzML")
            _copy_to_local(orig_mzml, local_mzml, log=_log)
            read_path = local_mzml

        # ---- build the DB from read_path (record the ORIGINAL mzML's identity) ----
        st = os.stat(orig_mzml)
        if os.path.exists(tmp):
            os.remove(tmp)
        conn = open_db(tmp)
        scan_meta: List[dict] = []
        n = 0
        reader = mzml.MzML(read_path, use_index=True)
        try:
            for spec in reader:
                spectrum = _extract_spectrum(spec)          # reuse mzml_utils
                extras = extract_extras(spec)
                insert_scan(conn, spectrum, extras)
                m = _HCD_RE.search(spectrum.filter_string or "")
                scan_meta.append(dict(
                    scan_num=spectrum.scan_num, activation_type=spectrum.activation_type,
                    iso_target=spectrum.isolation_window_target, precursor_mz=spectrum.precursor_mz,
                    sa=float(m.group(1)) if (m and spectrum.activation_type in ("EThcD", "ETD")) else None,
                ))
                n += 1
                if n % commit_every == 0:
                    conn.commit()
                    if progress:
                        progress(n, time.time() - t0)
        finally:
            try:
                reader.close()
            except Exception:
                pass
        n_pairs = finalize_pairs(conn, scan_meta)
        levels = {r["ms_level"]: r["c"] for r in conn.execute(
            "SELECT ms_level, COUNT(*) c FROM scans GROUP BY ms_level")}
        acts = {r["activation_type"]: r["c"] for r in conn.execute(
            "SELECT activation_type, COUNT(*) c FROM scans GROUP BY activation_type")}
        try:
            from . import __version__ as mu_ver
        except Exception:
            mu_ver = "unknown"
        set_meta(conn,
                 schema_version=SCHEMA_VERSION, source_mzml=orig_mzml,
                 source_size=st.st_size, source_mtime=int(st.st_mtime),
                 n_scans=n, n_pairs=n_pairs, ms_levels=levels, activations=acts,
                 build_seconds=round(time.time() - t0, 1), mzml_utils_version=mu_ver,
                 peak_dtype="float32", bit_exact=False)
        conn.commit()
        conn.close()
        conn = None
        os.replace(tmp, work)                # atomic; a crash never leaves a corrupt final DB
        if dst_on_net:
            _log(f"publishing DB to network: {db_path}")
            shutil.copy2(work, db_path)      # sequential copy of the finished DB to the network
        return dict(n_scans=n, n_pairs=n_pairs, ms_levels=levels, activations=acts,
                    seconds=round(time.time() - t0, 1), db_path=db_path,
                    db_bytes=os.path.getsize(db_path), local_copy=bool(local_mzml))
    finally:
        # remove ALL local temp files (source copy + local DB build file), success or failure.
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass
        if not keep_local_temp:
            leftovers = []
            if local_mzml:
                leftovers.append(local_mzml)
            if os.path.exists(tmp):          # leftover build temp on failure
                leftovers.append(tmp)
            if dst_on_net and os.path.exists(work) and os.path.abspath(work) != os.path.abspath(db_path):
                leftovers.append(work)       # local build copy of a network-dest DB
            for p in leftovers:
                try:
                    os.remove(p)
                    _log(f"removed local temp: {p}")
                except OSError:
                    pass


# --------------------------------------------------------------------------- #
# reader -- drop-in for mzml_utils.MzMLReader
# --------------------------------------------------------------------------- #
def _default_local_dir() -> str:
    return os.environ.get("SPECTRUM_CACHE_LOCAL", os.path.expanduser("~/.spectrum_cache_local"))


class SpectrumCache:
    """Reads the cache and returns mzml_utils.Spectrum objects.

    mirror_local: if the DB lives on a network share (/Volumes/...), copy it to a local
    dir once (small, one-time) and read from there -- avoids SMB random-read latency.
    in_memory: load the whole DB into RAM (one sequential read) and serve from memory.
    """

    def __init__(self, db_path: str, mirror_local: bool = True, in_memory: bool = False,
                 local_dir: Optional[str] = None, ms1_cap: Optional[int] = None):
        """``ms1_cap``: serve MS1 scans re-capped to their N most intense peaks (``cap_peaks``).
        On a cache patched from the RAW this reproduces the capped mzML exactly, so the same
        tool can be run on capped and uncapped MS1 side by side."""
        self.ms1_cap = ms1_cap
        self.source_path = db_path
        path = db_path
        if mirror_local and self._is_network(db_path):
            path = self._mirror(db_path, local_dir or _default_local_dir())
        self.path = path
        if in_memory:
            disk = sqlite3.connect(path)
            self._conn = sqlite3.connect(":memory:")
            disk.backup(self._conn)
            disk.close()
        else:
            self._conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        self._conn.row_factory = sqlite3.Row

    # -- network mirror ---------------------------------------------------- #
    @staticmethod
    def _is_network(p: str) -> bool:
        return _is_network(p)

    @staticmethod
    def _mirror(src: str, local_dir: str) -> str:
        os.makedirs(local_dir, exist_ok=True)
        dst = os.path.join(local_dir, os.path.basename(src))
        s = os.stat(src)
        if not (os.path.exists(dst) and os.stat(dst).st_size == s.st_size
                and int(os.stat(dst).st_mtime) >= int(s.st_mtime)):
            shutil.copy2(src, dst)  # one sequential read of the small DB
        return dst

    # -- MzMLReader-compatible API ---------------------------------------- #
    def _view(self, sp: Spectrum) -> Spectrum:
        if self.ms1_cap and sp.ms_level == 1 and len(sp.mz) > self.ms1_cap:
            sp.mz, sp.intensity = cap_peaks(sp.mz, sp.intensity, self.ms1_cap)
        return sp

    def get_spectrum(self, scan_num: int) -> Optional[Spectrum]:
        row = self._conn.execute("SELECT * FROM scans WHERE scan_num=?", (scan_num,)).fetchone()
        return self._view(_spectrum_from_row(row)) if row else None

    def iter_spectra(self, ms_level: Optional[int] = None,
                     rt_range: Optional[Tuple[float, float]] = None) -> Iterator[Spectrum]:
        """Iterate spectra in scan order, optionally only one MS level / RT window.

        The filters run in SQL, so a rejected scan is never zlib-decoded: an MS1-only
        sweep of a DDA run skips the MS2 scans that make up most of the file.
        `rt_range` is (min, max) in the stored unit (minutes). Same signature as
        MzMLReader.iter_spectra.
        """
        where, params = [], []
        if ms_level is not None:
            where.append("ms_level=?")
            params.append(int(ms_level))
        if rt_range is not None:
            where.append("rt BETWEEN ? AND ?")
            params.extend((float(rt_range[0]), float(rt_range[1])))
        sql = "SELECT * FROM scans"
        if where:
            sql += " WHERE " + " AND ".join(where)
        for row in self._conn.execute(sql + " ORDER BY scan_num", params):
            yield self._view(_spectrum_from_row(row))

    def find_best_ms1(self, ms2_scan: int, precursor_mz: float,
                      n_ms1: int = 3, tol_da: float = 0.02) -> Optional[int]:
        """MS1 scan with the strongest precursor peak among the preceding n_ms1 MS1s.

        Drop-in for MzMLReader.find_best_ms1 (same signature/semantics); uses the cache's
        fast get_spectrum lookups. Falls back to the closest preceding MS1 if the precursor
        peak is not detected.
        """
        ms1_scans: List[int] = []
        scan = ms2_scan - 1
        while len(ms1_scans) < n_ms1 and scan > 0:
            sp = self.get_spectrum(scan)
            if sp is not None and sp.ms_level == 1:
                ms1_scans.append(scan)
            scan -= 1
        best_scan, best_int = None, 0.0
        for ms1_scan in ms1_scans:
            sp = self.get_spectrum(ms1_scan)
            if sp is None:
                continue
            mask = np.abs(sp.mz - precursor_mz) < tol_da
            if mask.any():
                max_int = float(sp.intensity[mask].max())
                if max_int > best_int:
                    best_int, best_scan = max_int, ms1_scan
        if best_scan is None and ms1_scans:
            best_scan = ms1_scans[0]
        return best_scan

    def extract_metadata(self, scan_num: int) -> Optional[dict]:
        sp = self.get_spectrum(scan_num)
        if sp is None:
            return None
        return dict(scan_num=sp.scan_num, ms_level=sp.ms_level, rt=sp.rt,
                    filter_string=sp.filter_string, precursor_mz=sp.precursor_mz,
                    precursor_charge=sp.precursor_charge, activation_type=sp.activation_type,
                    isolation_window_target=sp.isolation_window_target,
                    isolation_window_lower=sp.isolation_window_lower,
                    isolation_window_upper=sp.isolation_window_upper,
                    isolation_width=sp.isolation_width,
                    ion_injection_time=sp.ion_injection_time)

    # -- extras + pairing -------------------------------------------------- #
    def raw_meta(self, scan_num: int) -> Optional[dict]:
        row = self._conn.execute("SELECT raw_meta_z FROM scans WHERE scan_num=?", (scan_num,)).fetchone()
        if not row or row["raw_meta_z"] is None:
            return None
        return json.loads(zlib.decompress(row["raw_meta_z"]).decode())

    def get_children(self, hcd_scan: int) -> List[int]:
        return [r["child_scan"] for r in self._conn.execute(
            "SELECT child_scan FROM pairs WHERE hcd_scan=? ORDER BY order_in_burst", (hcd_scan,))]

    def get_parent(self, child_scan: int) -> Optional[int]:
        r = self._conn.execute("SELECT hcd_scan FROM pairs WHERE child_scan=?", (child_scan,)).fetchone()
        return r["hcd_scan"] if r else None

    def find_preceding_ms1(self, scan_num: int) -> Optional[int]:
        """Master MS1 scan for a precursor -- O(1) via precursor_spectrum_ref, else the
        nearest earlier MS1."""
        r = self._conn.execute("SELECT precursor_spectrum_ref FROM scans WHERE scan_num=?",
                               (scan_num,)).fetchone()
        if r and r["precursor_spectrum_ref"]:
            return r["precursor_spectrum_ref"]
        r = self._conn.execute(
            "SELECT scan_num FROM scans WHERE ms_level=1 AND scan_num<? ORDER BY scan_num DESC LIMIT 1",
            (scan_num,)).fetchone()
        return r["scan_num"] if r else None

    def meta(self) -> dict:
        return {r["key"]: json.loads(r["value"]) for r in self._conn.execute("SELECT key, value FROM meta")}

    def is_stale(self, allow_network: bool = False) -> Optional[bool]:
        """True if the source mzML changed since build (size/mtime); None if unknown.

        SOURCE-BLIND by default: refuses to touch a network (/Volumes) source because
        os.stat / os.path.exists on a wedged SMB share can hang the process in
        uninterruptible I/O. Pass allow_network=True only for a source you know is
        reachable/local; otherwise this returns None (unknown) without touching it.
        """
        m = self.meta()
        src = m.get("source_mzml")
        if not src:
            return None
        if _is_network(src) and not allow_network:
            return None  # refuse to stat a network source -- can wedge on SMB
        if not os.path.exists(src):
            return None
        s = os.stat(src)
        return not (s.st_size == m.get("source_size") and int(s.st_mtime) == m.get("source_mtime"))

    def close(self):
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def _open_cache(db_path: str, in_memory=None, mirror_local=None, **kw):
    """Open a cache with defaults that keep local disk clean for NETWORK caches:
    a /Volumes/... DB loads into RAM (one sequential read, no local-disk copy); a local
    DB reads directly. Override with in_memory / mirror_local."""
    net = _is_network(db_path)
    if in_memory is None:
        in_memory = net
    if mirror_local is None:
        mirror_local = net and not in_memory
    return SpectrumCache(db_path, mirror_local=mirror_local, in_memory=in_memory, **kw)


def open_spectra(path: str, prefer_cache: bool = True, in_memory=None,
                 mirror_local=None, **kw):
    """Drop-in for MzMLReader(path). Given a .db/.sqlite/.cache path, opens the cache.
    Given an mzML path and prefer_cache=True, auto-uses the co-located cache
    (<mzml_dir>/spectra_cache/<stem>.spectra.db) when it exists, else MzMLReader.
    Network caches load into RAM by default so nothing is left on local disk.
    ``ms1_cap=N`` (caches only) serves MS1 re-capped to its N most intense peaks (``cap_peaks``),
    which on a cache patched from the RAW reproduces the capped mzML exactly."""
    p = str(path)
    if p.lower().endswith((".db", ".sqlite", ".cache")):
        return _open_cache(p, in_memory, mirror_local, **kw)
    if prefer_cache:
        cp = cache_path_for(p)
        if os.path.exists(cp):
            return _open_cache(cp, in_memory, mirror_local, **kw)
    from .reader import MzMLReader
    return MzMLReader(p)


# --------------------------------------------------------------------------- #
# MS1 from the RAW file: undo a per-scan peak cap (e.g. msconvert "1000 most intense")
# --------------------------------------------------------------------------- #
class PatchError(RuntimeError):
    """The RAW export does not match the cache; nothing was written."""


def _ms2_digest(conn: sqlite3.Connection) -> str:
    """sha256 over every non-MS1 row's scan number and peak blobs, in scan order."""
    import hashlib
    h = hashlib.sha256()
    for r in conn.execute("SELECT scan_num, mz_blob, int_blob FROM scans WHERE ms_level != 1 ORDER BY scan_num"):
        h.update(str(r[0]).encode())
        h.update(r[1] or b"")
        h.update(r[2] or b"")
    return h.hexdigest()


def cap_peaks(mz, intensity, n: int):
    """The ``n`` most intense peaks, excluding ties at the threshold (msconvert's rule), in m/z order.

    msconvert's "threshold count N most-intense" keeps the N most intense peaks and drops every
    peak tied at the threshold when keeping them all would exceed N, so a scan can end with N - 1.
    Applying this to a patched scan returns the capped scan exactly; ``patch_ms1`` checks that.
    """
    mz = np.asarray(mz)
    it = np.asarray(intensity)
    if n is None or len(it) <= n:
        return mz, it
    desc = np.sort(it)[::-1]
    thr = desc[n - 1]
    keep = it > thr if desc[n] == thr else it >= thr
    return mz[keep], it[keep]


def check_ms1_against_raw(cached_mz, cached_int, raw_mz, raw_int, mz_ppm: float = 0.5) -> Optional[str]:
    """None when the cached peaks are the RAW peaks after a most-intense cap; else the reason.

    The cached peaks must each match a distinct RAW peak (m/z within ``mz_ppm``, intensity equal
    at float32), and when the RAW scan has more peaks, the cached ones must be its most intense:
    no RAW peak left out may be more intense than the weakest one kept.
    """
    cmz = np.asarray(cached_mz, dtype=np.float64)
    cit = np.asarray(cached_int, dtype=np.float32)
    rmz = np.asarray(raw_mz, dtype=np.float64)
    rit = np.asarray(raw_int, dtype=np.float64)
    if len(cmz) == 0:
        return None if len(rmz) == 0 else f"cache scan is empty, the RAW scan has {len(rmz)} peaks"
    if len(rmz) == 0:
        return f"{len(cmz)} cached peaks, none in the RAW scan"
    if len(cmz) > len(rmz):
        return f"cache has more peaks ({len(cmz)}) than the RAW scan ({len(rmz)})"
    idx = np.clip(np.searchsorted(rmz, cmz), 1, len(rmz) - 1) if len(rmz) > 1 else np.zeros(len(cmz), int)
    if len(rmz) > 1:
        left, right = rmz[idx - 1], rmz[idx]
        j = np.where(np.abs(cmz - left) <= np.abs(cmz - right), idx - 1, idx)
    else:
        j = idx
    ppm = np.abs(rmz[j] - cmz) / cmz * 1e6
    if (ppm > mz_ppm).any():
        k = int(np.argmax(ppm))
        return f"cached peak m/z {cmz[k]:.5f} has no RAW peak within {mz_ppm} ppm (nearest {ppm[k]:.2f} ppm)"
    if len(np.unique(j)) != len(j):
        return "two cached peaks map to one RAW peak"
    rel = np.abs(rit[j].astype(np.float32) - cit) / np.maximum(np.abs(cit), 1.0)
    if (rel > 1e-6).any():
        k = int(np.argmax(rel))
        return f"intensity differs at m/z {cmz[k]:.5f}: cache {cit[k]:.6g}, RAW {rit[j][k]:.6g}"
    if len(j) < len(rmz):
        left_out = np.ones(len(rmz), bool)
        left_out[j] = False
        if rit[left_out].max() > float(cit.min()) * (1 + 1e-6):
            return (f"not a most-intense cap: a RAW peak left out ({rit[left_out].max():.6g}) is more "
                    f"intense than the weakest kept ({float(cit.min()):.6g})")
    return None


def merge_calibrated_ms1(cached_mz, cached_int, raw_mz, raw_int, max_shift_ppm: float = 50.0,
                         max_heldout_ppm: float = 1.0):
    """Uncap a scan whose kept peaks were m/z-recalibrated (e.g. MSFragger's *_calibrated.mzML).

    Each kept peak is paired with the RAW peak of identical intensity (float32) nearest in m/z,
    within ``max_shift_ppm``; every kept peak must pair, uniquely, and the kept peaks must be the
    RAW scan's most intense. The calibration shift (ppm) is then a function of m/z inside the scan
    and is interpolated from the pairs onto the peaks the cap removed (constant beyond the ends).
    Kept peaks keep their calibrated values exactly; added peaks get RAW intensity and
    interpolated calibrated m/z. Accuracy is measured on the scan itself, when peaks are added:
    shifts at odd-indexed pairs are predicted from the even-indexed ones, and the 99th percentile
    must stay under ``max_heldout_ppm``. A scan the cap did not cut is returned as cached.

    Returns (mz, intensity, info) or raises ValueError with the reason.
    """
    cmz = np.asarray(cached_mz, dtype=np.float64)
    cit = np.asarray(cached_int, dtype=np.float32)
    rmz = np.asarray(raw_mz, dtype=np.float64)
    rit = np.asarray(raw_int, dtype=np.float64)
    info = dict(n_kept=len(cmz), n_added=0, max_abs_shift_ppm=0.0, heldout_p99_ppm=0.0, heldout_max_ppm=0.0)
    if len(cmz) == 0:
        if len(rmz):
            raise ValueError(f"cache scan is empty, the RAW scan has {len(rmz)} peaks")
        return cmz, cit.astype(np.float64), info
    r32 = rit.astype(np.float32)
    by_int: Dict[float, List[int]] = {}
    for k, v in enumerate(r32):
        by_int.setdefault(float(v), []).append(k)
    used = set()
    pair = np.empty(len(cmz), dtype=int)
    for k, (m, v) in enumerate(zip(cmz, cit)):
        cand = [q for q in by_int.get(float(v), ()) if q not in used and abs(rmz[q] - m) / m * 1e6 <= max_shift_ppm]
        if not cand:
            raise ValueError(f"kept peak m/z {m:.5f} (intensity {v:.6g}) has no RAW peak of that intensity "
                             f"within {max_shift_ppm} ppm")
        q = min(cand, key=lambda q: abs(rmz[q] - m))
        used.add(q)
        pair[k] = q
    if len(pair) < len(rmz):
        left_out = np.ones(len(rmz), bool)
        left_out[pair] = False
        if rit[left_out].max() > float(cit.min()) * (1 + 1e-6):
            raise ValueError("not a most-intense cap: a RAW peak left out is more intense than the weakest kept")
    o = np.argsort(rmz[pair])
    a, b = rmz[pair][o], cmz[o]
    if (np.diff(b) < 0).any():
        raise ValueError("recalibration reorders peaks; not a smooth m/z function")
    shift = (b - a) / a * 1e6
    info["max_abs_shift_ppm"] = float(np.abs(shift).max())
    if len(pair) == len(rmz):          # nothing was cut, nothing to interpolate
        return cmz, cit.astype(np.float64), info
    if len(a) >= 10:
        tr_i, te_i = np.arange(0, len(a), 2), np.arange(1, len(a), 2)
        res = np.abs(np.interp(a[te_i], a[tr_i], shift[tr_i]) - shift[te_i])
        info["heldout_p99_ppm"] = float(np.percentile(res, 99))
        info["heldout_max_ppm"] = float(res.max())
        if info["heldout_p99_ppm"] > max_heldout_ppm:
            raise ValueError(f"calibration shift is not smooth in m/z (held-out p99 {info['heldout_p99_ppm']:.3f} ppm)")
    add = np.setdiff1d(np.arange(len(rmz)), pair)
    add_mz = rmz[add] * (1 + np.interp(rmz[add], a, shift) / 1e6)
    mz = np.concatenate([cmz, add_mz])
    it = np.concatenate([cit.astype(np.float64), rit[add]])
    srt = np.argsort(mz, kind="stable")
    info["n_added"] = int(len(add))
    return mz[srt], it[srt], info


def _default_raw_for(mzml_or_db: str) -> Optional[str]:
    """<mzml_dir>/<stem>.raw for an mzML (also stripping _calibrated / _uncalibrated)."""
    p = os.path.abspath(mzml_or_db)
    if p.lower().endswith(".db"):
        d = os.path.dirname(os.path.dirname(p))
        stem = os.path.basename(p)[: -len(".spectra.db")] if p.endswith(".spectra.db") else os.path.splitext(os.path.basename(p))[0]
    else:
        d, stem = os.path.dirname(p), os.path.splitext(os.path.basename(p))[0]
    for s in (stem, re.sub(r"_(un)?calibrated$", "", stem)):
        for ext in (".raw", ".RAW"):
            cand = os.path.join(d, s + ext)
            if os.path.isfile(cand):
                return cand
    return None


def patch_ms1(cache: str, raw_path: Optional[str] = None, *, export: Optional[dict] = None,
              pwiz_dir: Optional[str] = None, dry_run: bool = False, force: bool = False,
              mz_ppm: float = 0.5, rt_tol_min: float = 1e-4, mz_mode: str = "auto",
              log=None) -> dict:
    """Replace every MS1 scan's peaks in a spectrum cache with the full centroid list from the RAW.

    ``cache`` is the .spectra.db or the mzML it sits beside. ``raw_path`` defaults to the RAW next
    to the mzML. ``export`` ({scan: thermo_raw.RawMs1Scan}) bypasses the RAW read (tests).

    Checked before anything is written; any failure raises PatchError and leaves the cache as it
    was:
      * the cache and the RAW hold the same MS1 scan numbers, at the same retention times;
      * every cached MS1 peak is a RAW peak (m/z within ``mz_ppm``, same intensity), and where the
        RAW scan has more peaks the cached ones are its most intense, i.e. the only difference is
        a most-intense cap. A different centroider would fail here, so the patch never mixes two;
      * every RAW MS1 scan is centroid data (exporter source 1 or 2), not profile.
    ``mz_mode``: "exact" requires the cached m/z to be the RAW m/z and writes the RAW peaks.
    "calibrated" is for a cache built from a recalibrated mzML (MSFragger's *_calibrated.mzML):
    kept peaks stay as they are and the removed ones are added with m/z recalibrated by the
    scan's own shift curve (``merge_calibrated_ms1``). "auto" picks "calibrated" when the cache's
    source mzML name ends in _calibrated, else "exact".
    MS2 and higher rows are not touched: a digest of their peak blobs is compared before and after.
    TIC and base peak stay as recorded (the mzML takes them from the scan header).

    The cache records ``ms1_patch`` in its meta: RAW path, size and sha256, exporter version, the
    cap that was undone, peaks before and after, and the MS2 digest. A cache already patched from
    the same RAW is left alone; from a different RAW it raises unless ``force``.
    """
    from . import thermo_raw

    def _log(m):
        if log:
            log(m)

    db = cache if cache.lower().endswith((".db", ".sqlite", ".cache")) else cache_path_for(cache)
    if not os.path.isfile(db):
        raise FileNotFoundError(f"no spectrum cache at {db}; build it first")
    if export is None:
        raw_path = raw_path or _default_raw_for(cache)
        if not raw_path:
            raise FileNotFoundError(f"no RAW file found beside {cache}; pass raw_path / --raw")
    t0 = time.time()

    # Work on a local copy: SQLite writes over SMB are unsafe, and the original stays intact
    # until the patched copy is complete and checked.
    work_dir = tempfile.mkdtemp(prefix="ms1patch_", dir=_build_tmpdir())
    work = os.path.join(work_dir, os.path.basename(db))
    try:
        shutil.copy2(db, work)
        conn = sqlite3.connect(work)
        conn.row_factory = sqlite3.Row
        meta = {r["key"]: json.loads(r["value"]) for r in conn.execute("SELECT key, value FROM meta")}
        if mz_mode == "auto":
            src = os.path.splitext(os.path.basename(str(meta.get("source_mzml") or db)))[0]
            src = src[: -len(".spectra")] if src.endswith(".spectra") else src
            mz_mode = "calibrated" if src.endswith("_calibrated") else "exact"
        if mz_mode not in ("exact", "calibrated"):
            raise ValueError(f"mz_mode must be auto, exact or calibrated, not {mz_mode!r}")
        raw_id = None
        if raw_path:
            st = os.stat(raw_path)
            raw_id = dict(raw_path=os.path.abspath(raw_path), raw_size=st.st_size)
        prev = meta.get("ms1_patch")
        if prev:
            same = raw_id and prev.get("raw_size") == raw_id["raw_size"] and \
                os.path.basename(prev.get("raw_path", "")) == os.path.basename(raw_id["raw_path"])
            if same and not force:
                conn.close()
                _log(f"already patched from {prev.get('raw_path')}; nothing to do")
                return dict(prev, status="already_patched", db_path=db)
            if not force:
                conn.close()
                raise PatchError(f"{db} was patched from {prev.get('raw_path')}; pass force=True to redo it")

        if export is None:
            _log(f"exporting MS1 from {raw_path}")
            export = thermo_raw.export_ms1(raw_path, pwiz_dir=pwiz_dir)
        rows = conn.execute("SELECT scan_num, rt, n_peaks, mz_blob, int_blob FROM scans "
                            "WHERE ms_level = 1 ORDER BY scan_num").fetchall()
        cache_scans = [r["scan_num"] for r in rows]
        if set(cache_scans) != set(export):
            only_c = sorted(set(cache_scans) - set(export))[:5]
            only_r = sorted(set(export) - set(cache_scans))[:5]
            raise PatchError(f"MS1 scan numbers differ: cache-only {only_c}, RAW-only {only_r}")
        profile = [s for s, x in export.items() if x.source not in (1, 2)]
        if profile:
            raise PatchError(f"{len(profile)} RAW MS1 scans have no centroids (e.g. scan {profile[0]})")

        before_digest = _ms2_digest(conn)
        counts: Dict[int, int] = {}
        for r in rows:
            counts[r["n_peaks"]] = counts.get(r["n_peaks"], 0) + 1
        cap_hint = max(counts) if counts else None
        n_before = n_after = 0
        truncated, cap_values = 0, {}
        updates = []
        heldout = []
        max_shift = 0.0
        for r in rows:
            s = r["scan_num"]
            x = export[s]
            if r["rt"] is not None and abs(float(r["rt"]) - x.rt) > rt_tol_min:
                raise PatchError(f"scan {s}: retention time {r['rt']} in the cache, {x.rt} in the RAW")
            cmz, cit = decode_peaks(r["mz_blob"], r["int_blob"])
            if mz_mode == "exact":
                why = check_ms1_against_raw(cmz, cit, x.mz, x.intensity, mz_ppm=mz_ppm)
                if why:
                    raise PatchError(f"scan {s}: {why}" + (
                        " (if the cache was built from a recalibrated mzML, use mz_mode='calibrated')"
                        if "no RAW peak" in why else ""))
                new_mz, new_it = x.mz, x.intensity
            else:
                try:
                    new_mz, new_it, info = merge_calibrated_ms1(cmz, cit, x.mz, x.intensity)
                except ValueError as e:
                    raise PatchError(f"scan {s}: {e}")
                heldout.append(info["heldout_p99_ppm"])
                max_shift = max(max_shift, info["max_abs_shift_ppm"])
            if len(new_mz) > len(cmz):
                # Re-capping the patched scan must give back the cached scan, bit for bit at float32.
                rmz_c, rit_c = cap_peaks(np.asarray(new_mz, np.float32), np.asarray(new_it, np.float32), cap_hint)
                if len(rmz_c) != len(cmz) or not (np.array_equal(rmz_c, cmz.astype(np.float32))
                                                  and np.array_equal(rit_c, cit.astype(np.float32))):
                    raise PatchError(f"scan {s}: re-capping the patched scan at {cap_hint} does not "
                                     f"reproduce the cached scan ({len(rmz_c)} vs {len(cmz)} peaks)")
            n_before += len(cmz)
            n_after += len(new_mz)
            if len(x.mz) > len(cmz):
                truncated += 1
                cap_values[len(cmz)] = cap_values.get(len(cmz), 0) + 1
            mzb, itb = encode_peaks(new_mz, new_it)
            updates.append((mzb, itb, int(len(new_mz)), s))
        cap = max(cap_values, key=cap_values.get) if cap_values else None
        report = dict(
            db_path=db, raw_path=os.path.abspath(raw_path) if raw_path else None,
            raw_size=raw_id["raw_size"] if raw_id else None,
            exporter_version=thermo_raw.EXPORTER_VERSION, n_ms1=len(rows),
            ms1_truncated_before=truncated, cap_undone=cap, recap_reproduces_cache=True,
            truncated_peak_counts={str(k): v for k, v in sorted(cap_values.items())},
            ms1_peaks_before=n_before, ms1_peaks_after=n_after, ms2_digest=before_digest,
            mz_ppm=mz_ppm, rt_tol_min=rt_tol_min, mz_mode=mz_mode,
        )
        if mz_mode == "calibrated":
            report["max_abs_calibration_shift_ppm"] = round(max_shift, 3)
            report["heldout_p99_ppm_worst_scan"] = round(max(heldout), 4) if heldout else 0.0
            report["heldout_p99_ppm_median_scan"] = round(float(np.median(heldout)), 4) if heldout else 0.0
        if dry_run:
            conn.close()
            report["status"] = "checked (dry run, nothing written)"
            return report

        if raw_path:
            _log("hashing the RAW file")
            report["raw_sha256"] = thermo_raw.sha256_file(raw_path)
        conn.executemany("UPDATE scans SET mz_blob = ?, int_blob = ?, n_peaks = ? WHERE scan_num = ?", updates)
        if _ms2_digest(conn) != before_digest:
            raise PatchError("MS2 rows changed during the patch")
        report["patched_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        report["seconds"] = round(time.time() - t0, 1)
        set_meta(conn, ms1_patch={k: v for k, v in report.items() if k != "db_path"})
        conn.commit()
        conn.execute("VACUUM")
        conn.close()
        # Publish: copy beside the original, then an atomic rename over it.
        staged = db + ".ms1patch.tmp"
        shutil.copy2(work, staged)
        if os.path.getsize(staged) != os.path.getsize(work):
            os.remove(staged)
            raise PatchError("staged copy size mismatch; original left in place")
        os.replace(staged, db)
        report["status"] = "patched"
        report["db_bytes"] = os.path.getsize(db)
        _log(f"patched {len(rows)} MS1 scans ({truncated} were truncated at {cap}); "
             f"{n_before} -> {n_after} peaks")
        return report
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def ms1_cap_status(reader) -> dict:
    """Whether a reader's MS1 peaks look capped (a per-scan most-intense limit) or were patched.

    For a SpectrumCache: patched when its meta records ``ms1_patch``; otherwise capped when at
    least 5% of MS1 scans, and at least 3, share exactly the maximum peak count (a natural
    peak-count distribution has no such pile-up: its maximum is one scan). For an MzMLReader: read from the mzML header's dataProcessing
    ("most intense count"). Returns {capped, cap, patched, frac_at_cap, source}.
    """
    if isinstance(reader, SpectrumCache):
        m = reader.meta()
        if reader.ms1_cap:
            return dict(capped=True, cap=reader.ms1_cap, patched=bool(m.get("ms1_patch")),
                        frac_at_cap=float("nan"), source="re-capped view")
        if m.get("ms1_patch"):
            return dict(capped=False, cap=m["ms1_patch"].get("cap_undone"), patched=True,
                        frac_at_cap=0.0, source="cache-meta")
        tot = reader._conn.execute("SELECT COUNT(*) FROM scans WHERE ms_level = 1").fetchone()[0]
        top = reader._conn.execute("SELECT n_peaks, COUNT(*) FROM scans WHERE ms_level = 1 "
                                   "GROUP BY n_peaks ORDER BY n_peaks DESC LIMIT 1").fetchone()
        if not tot or not top:
            return dict(capped=False, cap=None, patched=False, frac_at_cap=0.0, source="peak-counts")
        frac = top[1] / tot
        capped = top[1] >= max(3, 0.05 * tot) and top[0] >= 50
        return dict(capped=capped, cap=int(top[0]) if capped else None, patched=False,
                    frac_at_cap=round(frac, 4), source="peak-counts")
    path = getattr(reader, "path", None) or getattr(reader, "filepath", None) or getattr(reader, "filename", None)
    cap = None
    if path and os.path.isfile(str(path)):
        with open(str(path), "rb") as fh:
            head = fh.read(200_000).decode("utf-8", "replace")
        m = re.search(r'most intense count[^"]*"\s+value="(\d+)"', head)
        if m:
            cap = int(m.group(1))
    return dict(capped=cap is not None, cap=cap, patched=False, frac_at_cap=float("nan"),
                source="mzml-header" if path else "unknown")


def warn_if_ms1_capped(reader, what: str = "MS1 peaks") -> dict:
    """Warn once per call site when MS1 is capped and not patched; returns ms1_cap_status."""
    import warnings
    st = ms1_cap_status(reader)
    if st["capped"]:
        warnings.warn(f"{what}: MS1 scans are capped at {st['cap']} peaks and the cache is not "
                      f"patched from the RAW (spectrum_cache patch-ms1); absent MS1 peaks may be "
                      f"below the cap, not missing.", stacklevel=2)
    return st


# --------------------------------------------------------------------------- #
# verify (run against the real mzML when the user chooses)
# --------------------------------------------------------------------------- #
def verify_cache(db_path: str, k: int = 25, seed_scans: Optional[List[int]] = None,
                 source: Optional[str] = None, allow_network: bool = False) -> dict:
    """Compare K scans between the cache and the source mzML via MzMLReader.

    Refuses a network (/Volumes) source by default -- reading a wedged SMB share can
    hang uninterruptibly. Pass an explicit LOCAL ``source`` override, or
    allow_network=True for a source you know is reachable. Raises RuntimeError rather
    than silently touching the share.
    """
    from .reader import MzMLReader
    cache = SpectrumCache(db_path, mirror_local=False)
    src = source or cache.meta().get("source_mzml")
    if src and _is_network(src) and not allow_network:
        cache.close()
        raise RuntimeError(
            f"refusing to verify against a network source ({src}); pass an explicit local "
            f"source= / --source or allow_network=True (reading a wedged SMB share can hang)")
    reader = MzMLReader(src)
    scans = seed_scans or [r["scan_num"] for r in cache._conn.execute(
        "SELECT scan_num FROM scans ORDER BY scan_num LIMIT ?", (k,))]
    mism = []
    for s in scans:
        a = cache.get_spectrum(s)
        b = reader.get_spectrum(s)
        if b is None:
            mism.append((s, "missing in mzML"))
            continue
        if a.n_peaks != b.n_peaks or (a.precursor_charge != b.precursor_charge):
            mism.append((s, f"n_peaks {a.n_peaks}/{b.n_peaks} charge {a.precursor_charge}/{b.precursor_charge}"))
            continue
        if a.n_peaks and not np.allclose(np.sort(a.mz), np.sort(b.mz), atol=1e-2):
            mism.append((s, "m/z differ"))
    reader.close()
    cache.close()
    return dict(checked=len(scans), mismatches=mism, ok=not mism)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _cli():
    ap = argparse.ArgumentParser(description="Spectrum annotation cache (build once, read fast).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build",
                       help="build a cache DB from an mzML (copy source LOCAL first, build LOCAL, copy DB to network)")
    b.add_argument("mzml")
    b.add_argument("-o", "--out", help="output .db path (may be on the network; build goes local then copies)")
    b.add_argument("--beside", action="store_true",
                   help="place cache in <mzml_dir>/spectra_cache/<stem>.spectra.db (built local, copied to network)")
    b.add_argument("--copy-to", help="also copy the finished DB here")
    b.add_argument("--direct", action="store_true",
                   help="read the mzML directly over the network (skip the default copy-to-local-first)")
    b.add_argument("--keep-local-temp", action="store_true",
                   help="keep the local temp files (source copy + build DB) instead of removing them")

    i = sub.add_parser("info", help="show cache provenance + counts (source-blind)")
    i.add_argument("db")
    i.add_argument("--check-source", action="store_true",
                   help="also stat the source mzML for staleness (do NOT use for a /Volumes source -- can hang)")

    v = sub.add_parser("verify", help="compare K scans against the source mzML")
    v.add_argument("db")
    v.add_argument("-k", type=int, default=25)
    v.add_argument("--source", help="verify against this LOCAL mzML instead of the recorded (maybe network) source")
    v.add_argument("--allow-network", action="store_true",
                   help="allow verifying against a /Volumes source (can hang on SMB)")

    pm = sub.add_parser("patch-ms1",
                        help="replace capped MS1 peaks in a cache with the full centroids from the RAW")
    pm.add_argument("target", nargs="+", help="spectrum cache .db or the mzML it sits beside")
    pm.add_argument("--raw", help="RAW file (default: <mzml_dir>/<stem>.raw); only with one target")
    pm.add_argument("--pwiz-dir", help="directory with the Thermo RawFileReader DLLs")
    pm.add_argument("--dry-run", action="store_true", help="check only, write nothing")
    pm.add_argument("--force", action="store_true", help="redo a cache patched from another RAW")
    pm.add_argument("--mz-ppm", type=float, default=0.5)
    pm.add_argument("--mz-mode", choices=["auto", "exact", "calibrated"], default="auto",
                    help="calibrated: the cache was built from a recalibrated mzML (default auto: by name)")

    ms = sub.add_parser("ms1-status", help="is MS1 capped, and was the cache patched?")
    ms.add_argument("target", nargs="+")

    args = ap.parse_args()
    if args.cmd == "patch-ms1":
        if args.raw and len(args.target) > 1:
            ap.error("--raw takes one target")
        bad = 0
        for tgt in args.target:
            try:
                res = patch_ms1(tgt, args.raw, pwiz_dir=args.pwiz_dir, dry_run=args.dry_run,
                                force=args.force, mz_ppm=args.mz_ppm, mz_mode=args.mz_mode,
                                log=lambda m: print(f"   {m}", flush=True))
                print(f"[{res['status']}] {tgt}: {res.get('n_ms1')} MS1, truncated "
                      f"{res.get('ms1_truncated_before')} at {res.get('cap_undone')}, peaks "
                      f"{res.get('ms1_peaks_before')} -> {res.get('ms1_peaks_after')}", flush=True)
            except (PatchError, FileNotFoundError, RuntimeError) as e:
                bad += 1
                print(f"[FAILED] {tgt}: {e}", flush=True)
        sys.exit(1 if bad else 0)
    if args.cmd == "ms1-status":
        for tgt in args.target:
            r = open_spectra(tgt)
            try:
                print(tgt, json.dumps(ms1_cap_status(r)))
            finally:
                r.close()
        return
    if args.cmd == "build":
        if args.beside:
            out = cache_path_for(args.mzml)
        else:
            out = args.out or (os.path.splitext(os.path.basename(args.mzml))[0] + ".spectra.db")
        mode = "direct network read" if args.direct else "copy-to-local-first (default)"
        print(f"[build] {args.mzml} -> {out}  [{mode}]")
        res = build_cache(args.mzml, out, keep_local_temp=args.keep_local_temp,
                          local_copy=(False if args.direct else None),
                          progress=lambda n, t: print(f"   {n} scans  {t:.0f}s", flush=True),
                          log=lambda m: print(f"   {m}", flush=True))
        print(f"[done] {res['n_scans']} scans, {res['n_pairs']} pairs, "
              f"{res['db_bytes']/1e6:.1f} MB, {res['seconds']}s"
              f"{' (via local copy)' if res.get('local_copy') else ''} -> {res['db_path']}")
        print(f"       ms_levels={res['ms_levels']} activations={res['activations']}")
        if args.copy_to:
            dst = args.copy_to
            if os.path.isdir(dst):
                dst = os.path.join(dst, os.path.basename(res["db_path"]))
            shutil.copy2(res["db_path"], dst)
            print(f"[copied] -> {dst}")
    elif args.cmd == "info":
        c = SpectrumCache(args.db, mirror_local=False)
        print(json.dumps(c.meta(), indent=2))
        if args.check_source:
            print("stale:", c.is_stale(allow_network=True))
        else:
            print("stale: (not checked -- source-blind; pass --check-source, but never for a /Volumes source)")
        c.close()
    elif args.cmd == "verify":
        res = verify_cache(args.db, k=args.k, source=args.source, allow_network=args.allow_network)
        print(json.dumps(res, indent=2, default=str))
        if not res.get("ok"):
            sys.exit(1)   # nonzero exit when the cache does not match the source


if __name__ == "__main__":
    _cli()
