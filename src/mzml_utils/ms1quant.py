"""Precursor MS1 intensity, IonQuant-like, for every PSM of a run -- and whether the cap censored it.

Ported from NGlycoMM ``data_analysis/tools/ms1_xic_extractor.py`` (2026-06-16), which was tuned
against IonQuant 1.11.20 ``Intensity`` on 4,887 precursors: Pearson(log) 0.924, Spearman 0.926,
median ratio 0.94, uniform across glycan classes. The rule, unchanged:

* For every MS1 scan within ``rttol`` minutes of the MS2 scan's retention time, sum the isotope
  cluster at the precursor m/z: isotopes ``iso_lo..iso_hi`` (relative to the anchor m/z) at
  spacing 1.0033548378 / z, each taken as the *nearest* peak within ``mztol`` ppm (on an exact
  tie, the higher-m/z peak). A scan counts when at least ``min_isotopes`` isotopes are found.
* A feature needs at least ``min_scans`` such scans; its intensity is the sum of the ``topn``
  largest scan sums ("top-3"). The apex scan sum, the trace area and the apex RT are reported too.

Parameters come from the search: ``Ms1QuantParams.from_workflow`` reads ``ionquant.mztol``,
``rttol``, ``minscans`` and ``minisotopes`` from a FragPipe ``fragpipe.workflow``.

Why not ``xic.extract_xics``: it makes one pass over the run per target and sums every peak in the
window. Here MS1 is indexed once per run and sliced by retention time, and each isotope takes the
nearest peak, the rule that was validated against IonQuant.

**Capped MS1.** On an mzML capped at N peaks per scan, a weak isotope can be missing because it fell
below the cap, not because it is absent. ``apex_floor_ratio`` is the strongest isotope of the apex
scan over that scan's lowest kept peak (inf when the scan is not at the cap). ``cap_censored`` is
True when that ratio is at most ``censor_factor`` (3), i.e. the feature itself sits near the cap
floor, or, for a missing feature, when any MS1 scan in its window sits at the cap. Weak isotopes
of a strong feature routinely sit near the floor and do not set the flag: they change the value
little, while a feature near the floor can appear, vanish or change scans when uncapped.
On a cache patched from the RAW (``spectrum_cache.patch_ms1``) nothing is capped and the flag is
always False. Run the same estimator on ``open_spectra(db, ms1_cap=N)`` to see the capped result.
"""
from __future__ import annotations

import configparser
import re
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Optional

import numpy as np

ISOTOPE_SPACING = 1.0033548378   # 13C - 12C


@dataclass(frozen=True)
class Ms1QuantParams:
    mztol_ppm: float = 10.0
    rttol_min: float = 0.4
    min_scans: int = 3
    min_isotopes: int = 2
    iso_lo: int = -1
    iso_hi: int = 4
    topn: int = 3
    censor_factor: float = 3.0

    @classmethod
    def from_workflow(cls, path: str, **overrides) -> "Ms1QuantParams":
        """IonQuant's mztol / rttol / minscans / minisotopes from a fragpipe.workflow."""
        text = open(path, encoding="utf-8", errors="replace").read()
        cp = configparser.ConfigParser(interpolation=None, strict=False)
        cp.read_string("[w]\n" + re.sub(r"(?m)^#.*$", "", text))
        w = cp["w"]
        kw = {}
        for key, name, typ in (("ionquant.mztol", "mztol_ppm", float), ("ionquant.rttol", "rttol_min", float),
                               ("ionquant.minscans", "min_scans", int), ("ionquant.minisotopes", "min_isotopes", int)):
            if key in w and w[key].strip():
                kw[name] = typ(float(w[key]))
        kw.update(overrides)
        return cls(**kw)


@dataclass
class Ms1Quant:
    top: float            # sum of the topn largest scan sums (0 when no feature)
    apex: float           # largest scan sum
    area: float           # trapezoid over the window's scan sums (scans without a cluster count 0)
    n_scans: int          # scans with >= min_isotopes isotopes
    apex_rt: float        # RT of the apex scan (nan when no feature)
    n_isotopes_apex: int  # isotopes found in the apex scan
    apex_floor_ratio: float  # strongest apex isotope / that scan's lowest kept peak (inf if not capped)
    cap_censored: bool    # see the module docstring

    def as_dict(self) -> dict:
        return asdict(self)


class Ms1Index:
    """Every MS1 scan of a run in memory, sorted by RT, plus every scan's RT.

    ``reader`` is anything from ``open_spectra`` (a cache, possibly re-capped with ``ms1_cap``,
    or an MzMLReader). The cap is read with ``spectrum_cache.ms1_cap_status``.
    """

    def __init__(self, reader):
        from .spectrum_cache import ms1_cap_status
        self.cap_status = ms1_cap_status(reader)
        cap = self.cap_status["cap"] if self.cap_status["capped"] else None
        self.scan_rt: Dict[int, float] = {}
        conn = getattr(reader, "_conn", None)
        if conn is not None:          # a SpectrumCache: every RT from SQL, MS1 peaks decoded only
            self.scan_rt = {int(s): float(rt) for s, rt in conn.execute("SELECT scan_num, rt FROM scans")}
            spectra = reader.iter_spectra(ms_level=1)
        else:
            spectra = self._ms1_and_rts(reader)
        rows = []
        for sp in spectra:
            mz = np.asarray(sp.mz, dtype=np.float64)
            it = np.asarray(sp.intensity, dtype=np.float64)
            if len(mz) and np.any(np.diff(mz) < 0):
                o = np.argsort(mz, kind="stable")
                mz, it = mz[o], it[o]
            at_cap = bool(cap) and len(mz) >= cap - 1      # N-1 when ties at the threshold were dropped
            floor = float(it.min()) if (at_cap and len(it)) else 0.0
            rows.append((float(sp.rt), mz, it, floor))
        # RT in seconds (some converters) -> minutes, as the NGlycoMM extractor did.
        rt_max = max([r[0] for r in rows] + list(self.scan_rt.values()) + [0.0])
        if rt_max > 600:
            rows = [(r[0] / 60.0,) + r[1:] for r in rows]
            self.scan_rt = {k: v / 60.0 for k, v in self.scan_rt.items()}
        rows.sort(key=lambda r: r[0])
        self.rt = np.array([r[0] for r in rows], dtype=np.float64)
        self.mz = [r[1] for r in rows]
        self.intensity = [r[2] for r in rows]
        self.floor = np.array([r[3] for r in rows], dtype=np.float64)
        self.capped = cap is not None

    def _ms1_and_rts(self, reader):
        for sp in reader.iter_spectra():
            self.scan_rt[int(sp.scan_num)] = float(sp.rt)
            if sp.ms_level == 1:
                yield sp

    def __len__(self):
        return len(self.rt)


def _nearest_peak(mz: np.ndarray, target: float, tol_da: float) -> int:
    """Index of the peak nearest ``target`` within ``tol_da`` (ties: the higher m/z), else -1."""
    k = int(np.searchsorted(mz, target))
    best, bestd = -1, tol_da
    for kk in (k - 1, k):
        if 0 <= kk < mz.shape[0]:
            d = abs(float(mz[kk]) - target)
            if d <= bestd:
                bestd, best = d, kk
    return best


def isotope_cluster(mz: np.ndarray, intensity: np.ndarray, anchor_mz: float, charge: int,
                    params: Ms1QuantParams = Ms1QuantParams()):
    """(sum, n_found, strongest_found) of the isotope cluster in one scan."""
    s, n, strongest = 0.0, 0, 0.0
    for i in range(params.iso_lo, params.iso_hi + 1):
        target = anchor_mz + i * ISOTOPE_SPACING / charge
        j = _nearest_peak(mz, target, target * params.mztol_ppm * 1e-6)
        if j >= 0:
            v = float(intensity[j])
            s += v
            n += 1
            strongest = max(strongest, v)
    return s, n, strongest


def precursor_ms1(index: Ms1Index, anchor_mz: float, charge: int, rt_center: float,
                  params: Ms1QuantParams = Ms1QuantParams()) -> Ms1Quant:
    """The feature at ``anchor_mz`` / ``charge`` around ``rt_center`` (minutes)."""
    nan = float("nan")
    if not charge or charge < 1 or not (anchor_mz > 0) or not (rt_center == rt_center):
        return Ms1Quant(0.0, 0.0, 0.0, 0, nan, 0, nan, False)
    i0 = int(np.searchsorted(index.rt, rt_center - params.rttol_min))
    i1 = int(np.searchsorted(index.rt, rt_center + params.rttol_min))
    sums, rts, kept = [], [], []
    window_at_cap = False
    for idx in range(i0, i1):
        s, n, strongest = isotope_cluster(index.mz[idx], index.intensity[idx], anchor_mz, charge, params)
        window_at_cap |= index.floor[idx] > 0
        ok = n >= params.min_isotopes and s > 0
        rts.append(index.rt[idx])
        sums.append(s if ok else 0.0)
        if ok:
            kept.append((s, idx, n, strongest))
    if len(kept) < params.min_scans:
        return Ms1Quant(0.0, 0.0, 0.0, len(kept), nan, 0, nan, bool(window_at_cap))
    vals = np.array([k[0] for k in kept])
    top = float(np.sort(vals)[-params.topn:].sum())
    a = int(np.argmax(vals))
    _, aidx, an, astrong = kept[a]
    floor = index.floor[aidx]
    ratio = float(astrong / floor) if floor > 0 else float("inf")
    trapz = getattr(np, "trapezoid", None) or np.trapz
    area = float(trapz(sums, rts)) if len(rts) > 1 else 0.0
    return Ms1Quant(top, float(vals[a]), area, len(kept), float(index.rt[aidx]), int(an), ratio,
                    bool(ratio <= params.censor_factor))


def quantify_precursors(source, rows: Iterable[dict], params: Ms1QuantParams = Ms1QuantParams(),
                        ms1_cap: Optional[int] = None) -> List[dict]:
    """MS1 for many precursors of one run.

    ``source`` is a path (opened with ``open_spectra``, honouring ``ms1_cap``) or an open reader.
    Each row needs ``ms2_scan`` and ``charge``; ``anchor_mz`` defaults to the MS2 scan's selected
    ion m/z from the reader (instrument coordinates) and ``rt_center`` to the MS2 scan's RT.
    Returns one dict per row: the row's keys plus the ``Ms1Quant`` fields and ``anchor_source``.
    """
    from .spectrum_cache import open_spectra
    own = not hasattr(source, "iter_spectra")
    reader = open_spectra(str(source), **({"ms1_cap": ms1_cap} if ms1_cap else {})) if own else source
    try:
        index = Ms1Index(reader)
        out = []
        for r in rows:
            r = dict(r)
            scan = int(r["ms2_scan"])
            anchor = r.get("anchor_mz")
            src = "given"
            if anchor is None or not (anchor == anchor):
                sp = reader.get_spectrum(scan)
                anchor = float(sp.precursor_mz) if sp is not None else float("nan")
                src = "selected_ion"
            rtc = r.get("rt_center")
            if rtc is None or not (rtc == rtc):
                rtc = index.scan_rt.get(scan, float("nan"))
            q = precursor_ms1(index, float(anchor), int(r["charge"]), float(rtc), params)
            r.update(q.as_dict(), anchor_mz=float(anchor), anchor_source=src, rt_center=float(rtc))
            out.append(r)
        return out
    finally:
        if own:
            reader.close()
