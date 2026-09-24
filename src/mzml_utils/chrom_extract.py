"""
Table-driven chromatogram extraction across many runs.

The Python counterpart of ``chromExtract()`` in the Bioconductor package
Chromatograms (Louail, Gatto, Gibb & Rainer, Anal Chem 2026). A peak table from any
software -- MSFragger PSMs, a published supplementary table, a hand-written target
list -- names an m/z window and a retention-time window per row. Each row's
chromatogram is pulled from the raw data, and every column of the row is carried
through unchanged, so the result stays traceable to the table that asked for it.

What differs from :func:`~mzml_utils.xic.extract_xics`:

* every row has its OWN retention-time window (and m/z window), so glycoforms that
  elute minutes apart are extracted in one pass without a whole-run trace each;
* one pass per run serves all rows of that run: each MS1 scan is decoded once, and
  all active rows are summed in one vectorised ``searchsorted`` over a cumulative
  intensity sum instead of a mask per target;
* many runs in one call, with ``by=`` restricting a row to its own run (a PSM
  belongs to the file it was identified in) or, without ``by``, every row extracted
  from every run (cross-sample overlays: one glycoform across all fractions);
* optional :func:`~mzml_utils.xic.peak_metrics` per chromatogram.

Typical use::

    from mzml_utils import extract_chromatograms, chromatogram_rows

    table = [
        {"psg": "P01024_N1115_N4H5F1A1", "mz": 1204.564, "rt_min": 61.5, "rt_max": 63.5,
         "run": "fraction_07", "ms2_rt": 62.4},
        ...
    ]
    res = extract_chromatograms(table, ["/data/fraction_07.mzML", ...], by="run",
                                apex_rt="ms2_rt", nearest_valley=True)
    rows = chromatogram_rows(res)      # list of flat dicts; pandas.DataFrame(rows)
"""

from __future__ import annotations

import inspect
import os
import warnings
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .xic import XIC, PeakMetrics, peak_metrics

SourceSpec = Union[Sequence[Union[str, os.PathLike]], Mapping[str, Any]]

_SOURCE_SUFFIXES = (".spectra.db", ".mzml", ".mzxml", ".db")


def source_key(path: Union[str, os.PathLike]) -> str:
    """The run name of a spectrum source: its file name without the mzML/cache suffix.

    ``/data/spectra_cache/run_01.spectra.db`` and ``/data/run_01.mzML`` both give
    ``run_01``; this is what a ``by=`` column must contain.
    """
    name = os.path.basename(str(path))
    low = name.lower()
    for suf in _SOURCE_SUFFIXES:
        if low.endswith(suf):
            return name[: len(name) - len(suf)]
    return name


@dataclass
class ExtractedChromatogram:
    """One row of the peak table, extracted from one run."""

    row: Dict[str, Any]
    """The peak-table row, carried through unchanged."""
    source: str
    """Run name (:func:`source_key`, or the key given in a ``sources`` mapping)."""
    xic: XIC
    metrics: Optional[PeakMetrics] = None


def _normalize_sources(sources: SourceSpec) -> List[Tuple[str, Any]]:
    if isinstance(sources, Mapping):
        return [(str(k), v) for k, v in sources.items()]
    if isinstance(sources, (str, os.PathLike)):
        sources = [sources]
    out = [(source_key(p), p) for p in sources]
    keys = [k for k, _ in out]
    if len(set(keys)) != len(keys):
        dup = sorted({k for k in keys if keys.count(k) > 1})
        raise ValueError(f"sources share a run name {dup}; pass a {{name: path}} mapping")
    return out


def _row_windows(rows: List[Mapping[str, Any]], mz: str, mz_min: str, mz_max: str,
                 rt_min: str, rt_max: str, tolerance: float, unit: str):
    """Arrays (mz_lo, mz_hi, rt_lo, rt_hi, target_mz), validated like chromExtract."""
    if unit not in ("ppm", "Da"):
        raise ValueError(f"Unknown unit '{unit}'. Use 'ppm' or 'Da'.")
    n = len(rows)
    lo = np.empty(n)
    hi = np.empty(n)
    target = np.empty(n)
    rlo = np.empty(n)
    rhi = np.empty(n)
    for i, r in enumerate(rows):
        try:
            rlo[i] = float(r[rt_min])
            rhi[i] = float(r[rt_max])
        except KeyError as e:
            raise ValueError(f"peak table row {i} has no {e.args[0]!r} column") from None
        if np.isnan(rlo[i]) or np.isnan(rhi[i]):
            raise ValueError(f"peak table row {i}: {rt_min}/{rt_max} is missing (NaN)")
        if rlo[i] > rhi[i]:
            raise ValueError(f"peak table row {i}: {rt_min} > {rt_max}")
        has_min, has_max = mz_min in r, mz_max in r
        if has_min != has_max:
            raise ValueError(f"peak table row {i}: give both {mz_min} and {mz_max}, or neither")
        if has_min:
            lo[i], hi[i] = float(r[mz_min]), float(r[mz_max])
            target[i] = float(r[mz]) if mz in r else (lo[i] + hi[i]) / 2.0
        else:
            if mz not in r:
                raise ValueError(f"peak table row {i} has neither {mz!r} nor {mz_min!r}/{mz_max!r}")
            t = float(r[mz])
            tol = t * tolerance / 1e6 if unit == "ppm" else tolerance
            lo[i], hi[i], target[i] = t - tol, t + tol, t
        if np.isnan(lo[i]) or np.isnan(hi[i]):
            raise ValueError(f"peak table row {i}: m/z is missing (NaN)")
    return lo, hi, rlo, rhi, target


def _iter_pushed(reader, ms_level: Optional[int], rt_range: Tuple[float, float]):
    """iter_spectra with the MS-level / RT filters pushed down when accepted."""
    accepted = inspect.signature(reader.iter_spectra).parameters
    pushed = {k: v for k, v in (("ms_level", ms_level), ("rt_range", rt_range))
              if v is not None and k in accepted}
    for spec in reader.iter_spectra(**pushed):
        if ms_level is not None and spec.ms_level != ms_level:
            continue
        if not (rt_range[0] <= spec.rt <= rt_range[1]):
            continue
        yield spec


def _extract_one_source(reader, idx: np.ndarray, lo, hi, rlo, rhi, ms_level):
    """(rt, intensity) per row index in ``idx`` from one open reader, in one pass."""
    s_lo, s_hi, s_rlo, s_rhi = lo[idx], hi[idx], rlo[idx], rhi[idx]
    rt_parts: List[np.ndarray] = []
    row_parts: List[np.ndarray] = []
    val_parts: List[np.ndarray] = []
    for spec in _iter_pushed(reader, ms_level, (float(s_rlo.min()), float(s_rhi.max()))):
        rt = float(spec.rt)
        active = np.flatnonzero((s_rlo <= rt) & (rt <= s_rhi))
        if active.size == 0:
            continue
        mz_a = np.asarray(spec.mz, dtype=float)
        int_a = np.asarray(spec.intensity, dtype=float)
        if mz_a.size and np.any(np.diff(mz_a) < 0):
            order = np.argsort(mz_a, kind="stable")
            mz_a, int_a = mz_a[order], int_a[order]
        csum = np.concatenate(([0.0], np.cumsum(int_a)))
        a = np.searchsorted(mz_a, s_lo[active], side="left")
        b = np.searchsorted(mz_a, s_hi[active], side="right")
        val_parts.append(csum[b] - csum[a])
        row_parts.append(active)
        rt_parts.append(np.full(active.size, rt))

    out: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    if row_parts:
        rows = np.concatenate(row_parts)
        rts = np.concatenate(rt_parts)
        vals = np.concatenate(val_parts)
        order = np.lexsort((rts, rows))  # by row, then retention time (stable)
        rows, rts, vals = rows[order], rts[order], vals[order]
        cuts = np.flatnonzero(np.diff(rows)) + 1
        for r_rows, r_rt, r_val in zip(np.split(rows, cuts), np.split(rts, cuts),
                                       np.split(vals, cuts)):
            out[int(idx[r_rows[0]])] = (r_rt, r_val)
    return out


def extract_chromatograms(peak_table: Iterable[Mapping[str, Any]],
                          sources: SourceSpec, *,
                          by: Optional[str] = None,
                          mz: str = "mz",
                          mz_min: str = "mz_min",
                          mz_max: str = "mz_max",
                          rt_min: str = "rt_min",
                          rt_max: str = "rt_max",
                          tolerance: float = 10.0,
                          unit: str = "ppm",
                          ms_level: Optional[int] = 1,
                          metrics: bool = True,
                          apex_rt: Optional[str] = None,
                          **boundary_kwargs) -> List[ExtractedChromatogram]:
    """Extract one chromatogram per (peak-table row, run), carrying the row through.

    Args:
        peak_table: Rows as mappings (``df.to_dict("records")`` for a DataFrame). Each
            row needs ``rt_min``/``rt_max`` and either ``mz`` (window = +/- ``tolerance``)
            or both ``mz_min``/``mz_max``. A missing or NaN window is an error, as in
            chromExtract -- a row is never silently extracted over the whole run.
        sources: Spectrum files (mzML or ``.spectra.db``; opened with
            :func:`~mzml_utils.open_spectra`, so a co-located cache is used), or a
            ``{run_name: path_or_open_reader}`` mapping. Readers passed in are not closed.
        by: Row column holding the run name. When given, a row is extracted only from
            the run whose name equals it; a run name matching no source is an error.
            When ``None``, every row is extracted from every source.
        mz, mz_min, mz_max, rt_min, rt_max: Column names.
        tolerance, unit: m/z window for rows given by ``mz`` (``'ppm'`` or ``'Da'``).
        ms_level: MS level extracted (1 = precursor traces). ``None`` uses every scan.
        metrics: Compute :func:`~mzml_utils.xic.peak_metrics` for each chromatogram.
        apex_rt: Optional column of retention times the measured peak must contain
            (e.g. the PSM's MS2 scan time); passed as ``contains_rt`` to
            :func:`~mzml_utils.xic.peak_boundary`. Without it the tallest peak in the
            window is measured. A row whose value is missing falls back to the tallest.
        **boundary_kwargs: ``threshold``, ``baseline_threshold``, ``baseline_quantile``,
            ``nearest_valley`` passed to :func:`~mzml_utils.xic.peak_boundary`.
            ``nearest_valley=True`` is recommended for these zero-filled traces.

    Returns:
        One :class:`ExtractedChromatogram` per (row, run) pair, in row order then source
        order. A row whose window holds no scan gets an empty XIC and empty metrics.
        Intensities within an m/z window are summed (the ``mode='sum'`` of
        ``extract_xics``); they agree with ``extract_xics`` to floating-point rounding.
    """
    rows = [dict(r) for r in peak_table]
    srcs = _normalize_sources(sources)
    if not rows or not srcs:
        return []
    lo, hi, rlo, rhi, target = _row_windows(rows, mz, mz_min, mz_max, rt_min, rt_max,
                                            tolerance, unit)

    names = [k for k, _ in srcs]
    if by is not None:
        missing = [i for i, r in enumerate(rows) if by not in r]
        if missing:
            raise ValueError(f"peak table rows {missing[:5]} have no {by!r} column")
        unknown = sorted({str(r[by]) for r in rows} - set(names))
        if unknown:
            raise ValueError(f"{by!r} values match no source: {unknown[:5]} "
                             f"(sources: {names[:5]}{' ...' if len(names) > 5 else ''})")

    traces: Dict[Tuple[int, str], Tuple[np.ndarray, np.ndarray]] = {}
    for name, src in srcs:
        if by is None:
            idx = np.arange(len(rows))
        else:
            idx = np.array([i for i, r in enumerate(rows) if str(r[by]) == name], dtype=int)
        if idx.size == 0:
            continue
        if hasattr(src, "iter_spectra"):
            got = _extract_one_source(src, idx, lo, hi, rlo, rhi, ms_level)
        else:
            from .spectrum_cache import open_spectra
            with open_spectra(str(src)) as reader:
                got = _extract_one_source(reader, idx, lo, hi, rlo, rhi, ms_level)
        for i, tr in got.items():
            traces[(i, name)] = tr

    empty = (np.array([], dtype=float), np.array([], dtype=float))
    out: List[ExtractedChromatogram] = []
    for i, row in enumerate(rows):
        for name in names:
            if by is not None and str(row[by]) != name:
                continue
            rt_arr, int_arr = traces.get((i, name), empty)
            xic = XIC(target_mz=float(target[i]), rt=rt_arr, intensity=int_arr,
                      name=f"{name}:{i}", ms_level=ms_level if ms_level is not None else 0)
            pm = None
            if metrics:
                t = row.get(apex_rt) if apex_rt is not None else None
                t = None if t is None or np.isnan(float(t)) else float(t)
                pm = peak_metrics(xic, contains_rt=t, **boundary_kwargs)
            out.append(ExtractedChromatogram(row=row, source=name, xic=xic, metrics=pm))
    return out


def chromatogram_rows(results: Sequence[ExtractedChromatogram]) -> List[Dict[str, Any]]:
    """Flatten results to one dict per chromatogram: the row's own columns, then
    ``source``, ``n_scans`` and the metric columns (``pandas.DataFrame(rows)`` ready).

    A row column that shares a name with an added column is replaced, with a warning
    (chromExtract does the same).
    """
    flat: List[Dict[str, Any]] = []
    warned = set()
    for res in results:
        d = dict(res.row)
        added: Dict[str, Any] = {"source": res.source, "n_scans": res.xic.n_points}
        if res.metrics is not None:
            added.update(res.metrics.as_dict())
        for k, v in added.items():
            if k in d and k not in warned:
                warnings.warn(f"chromatogram_rows: peak-table column {k!r} replaced by the "
                              f"computed value", stacklevel=2)
                warned.add(k)
            d[k] = v
        flat.append(d)
    return flat
