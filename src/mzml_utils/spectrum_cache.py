#!/usr/bin/env python3
"""
spectrum-cache: build-once, read-fast local cache of everything spectrum annotation needs.

Reads an mzML ONE time (the only slow / network step) and writes a compact, lossless,
indexed SQLite cache. All later annotation reads go through ``SpectrumCache``, which
returns ``mzml_utils.Spectrum`` objects -- a drop-in for ``mzml_utils.MzMLReader`` -- so
both the ``nglyco-analysis`` and ``oglyco-analysis`` skills work unchanged via
``open_spectra()``.

Design (see PLAN):
  * LOSSLESS raw peaks (float32 + zlib). N-glyco matches raw (un-deisotoped) spectra and
    O-glyco deisotopes on the fly -- a filtered/deisotoped cache would break N-glyco.
  * Captures MORE than the 15-field Spectrum dataclass: ion injection time, FAIMS CV,
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
    mz_blob               BLOB,        -- zlib(float32[])  RAW, lossless
    int_blob              BLOB,        -- zlib(float32[])  RAW, lossless
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
    """Pull the mzML fields the 15-field Spectrum dataclass does NOT carry."""
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

    return (
        int(spectrum.scan_num), int(spectrum.ms_level), spectrum.activation_type,
        num(spectrum.rt), spectrum.filter_string,
        num(spectrum.precursor_mz), int(spectrum.precursor_charge or 0),
        num(spectrum.precursor_intensity), extras.get("precursor_spectrum_ref"),
        num(spectrum.isolation_window_target), num(spectrum.isolation_window_lower),
        num(spectrum.isolation_window_upper), num(extras.get("ion_injection_time")),
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
                 lossless=True)
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
                 local_dir: Optional[str] = None):
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
    def get_spectrum(self, scan_num: int) -> Optional[Spectrum]:
        row = self._conn.execute("SELECT * FROM scans WHERE scan_num=?", (scan_num,)).fetchone()
        return _spectrum_from_row(row) if row else None

    def iter_spectra(self) -> Iterator[Spectrum]:
        for row in self._conn.execute("SELECT * FROM scans ORDER BY scan_num"):
            yield _spectrum_from_row(row)

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
                    isolation_width=sp.isolation_width)

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

    def is_stale(self) -> Optional[bool]:
        """True if the source mzML changed since build (size/mtime); None if unknown."""
        m = self.meta()
        src = m.get("source_mzml")
        if not src or not os.path.exists(src):
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
    Network caches load into RAM by default so nothing is left on local disk."""
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
# verify (run against the real mzML when the user chooses)
# --------------------------------------------------------------------------- #
def verify_cache(db_path: str, k: int = 25, seed_scans: Optional[List[int]] = None) -> dict:
    """Compare K scans between the cache and the source mzML via MzMLReader."""
    from .reader import MzMLReader
    cache = SpectrumCache(db_path, mirror_local=False)
    src = cache.meta().get("source_mzml")
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

    i = sub.add_parser("info", help="show cache provenance + counts")
    i.add_argument("db")

    v = sub.add_parser("verify", help="compare K scans against the source mzML")
    v.add_argument("db")
    v.add_argument("-k", type=int, default=25)

    args = ap.parse_args()
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
        print("stale:", c.is_stale())
        c.close()
    elif args.cmd == "verify":
        print(json.dumps(verify_cache(args.db, k=args.k), indent=2, default=str))


if __name__ == "__main__":
    _cli()
