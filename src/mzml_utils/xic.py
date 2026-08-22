"""
Extracted ion chromatogram (XIC) extraction.

An XIC is a trace of ion intensity versus retention time for one target m/z.
This module builds XICs from any spectrum source that yields ``Spectrum``
objects (an :class:`~mzml_utils.reader.MzMLReader`, a
:class:`~mzml_utils.spectrum_cache.SpectrumCache`, or a path handed to
``open_spectra``) in a single pass over the file.

Two extraction domains are covered by the ``ms_level`` argument:

* ``ms_level=1`` -- precursor XICs. Extract a target precursor/(glyco)peptide
  m/z from every MS1 scan (e.g. tracking co-eluting glycoforms).
* ``ms_level=2`` -- diagnostic/fragment XICs. Extract a target fragment m/z
  (e.g. oxonium ions 138.055, 204.087, 274.092) from every MS2 scan, giving
  an overview of where glycopeptides elute. Optionally restrict to MS2 scans
  of a single precursor (``precursor_mz``) or a single activation (``activation``).

Retention time is reported in the reader's native unit (minutes for Thermo
msconvert-derived mzML, matching the rest of this ecosystem). ``rt_range`` is
interpreted in the same unit.

Typical use::

    from mzml_utils import extract_xics, OXONIUM_IONS

    # Oxonium-ion overview across the whole run (Fig S1 top-block style)
    cs = extract_xics("sample.mzML",
                      {k: OXONIUM_IONS[k] for k in ("HexNAc", "NeuAc")},
                      ms_level=2, tolerance=20.0, unit="ppm", activation="HCD")
    print(cs["HexNAc"].max_intensity, cs["HexNAc"].apex_rt)

    # MS1 precursor XIC for one glycoform over a retention window
    from mzml_utils import extract_xic
    xic = extract_xic("sample.mzML", 1204.564, ms_level=1, rt_range=(62.5, 74.1))
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

# numpy>=2 renamed trapz -> trapezoid; keep working on both.
# NB: must be a conditional, not getattr(np, "trapezoid", getattr(np, "trapz")) -- the default
# argument there is evaluated eagerly, so on numpy>=2 (where trapz is gone) it raises
# AttributeError before the trapezoid result can ever be used. That inverted the intent: it
# worked on numpy 1 and broke on numpy 2.
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

TargetSpec = Union[Mapping[str, float], Sequence[Union[float, Tuple[str, float]]]]


@dataclass
class XIC:
    """A single extracted ion chromatogram (intensity versus retention time)."""

    target_mz: float
    rt: np.ndarray
    intensity: np.ndarray
    name: str = ""
    ms_level: int = 0

    @property
    def n_points(self) -> int:
        return int(len(self.rt))

    @property
    def max_intensity(self) -> float:
        """Peak intensity of the trace -- the Thermo 'NL' (normalized level) value."""
        return float(self.intensity.max()) if len(self.intensity) else 0.0

    @property
    def apex_rt(self) -> float:
        """Retention time of the most intense point."""
        if not len(self.intensity):
            return 0.0
        return float(self.rt[int(np.argmax(self.intensity))])

    @property
    def area(self) -> float:
        """Trapezoidal peak area over retention time."""
        if len(self.rt) < 2:
            return 0.0
        return float(_trapz(self.intensity, self.rt))

    def _half_max_span(self):
        """(apex_idx, half, lo_idx, hi_idx) of the CONTIGUOUS half-maximum region.

        Contiguous on purpose: an XIC often contains more than one peak (a second
        elution, or an interfering species at the same m/z), and counting every point
        above half maximum anywhere in the trace would merge them into one absurdly
        wide peak. Only the region containing the apex is the peak.
        """
        n = len(self.intensity)
        if n == 0:
            return None
        apex = int(np.argmax(self.intensity))
        half = float(self.intensity[apex]) / 2.0
        if half <= 0.0:
            return None
        lo = apex
        while lo > 0 and self.intensity[lo - 1] >= half:
            lo -= 1
        hi = apex
        while hi < n - 1 and self.intensity[hi + 1] >= half:
            hi += 1
        return apex, half, lo, hi

    @property
    def fwhm(self) -> float:
        """Full width at half maximum, in the retention-time unit of `rt` (minutes).

        The crossings are linearly interpolated between the bracketing scans rather
        than snapped to them, because at 9-15 points across a peak, snapping
        quantises the width in steps of a whole cycle time.

        Returns 0.0 when there is no determinable peak (empty trace, flat, or a
        single point). Check `fwhm_is_truncated` before treating the value as a
        measurement: a peak running off either end of the trace yields a LOWER BOUND.
        """
        span = self._half_max_span()
        if span is None or len(self.rt) < 2:
            return 0.0
        _apex, half, lo, hi = span

        def cross(i_in, i_out):
            """RT where the trace crosses `half` between an inside and outside point."""
            y_in, y_out = float(self.intensity[i_in]), float(self.intensity[i_out])
            x_in, x_out = float(self.rt[i_in]), float(self.rt[i_out])
            if y_in == y_out:
                return x_out
            return x_out + (half - y_out) * (x_in - x_out) / (y_in - y_out)

        left = cross(lo, lo - 1) if lo > 0 else float(self.rt[lo])
        right = cross(hi, hi + 1) if hi < len(self.rt) - 1 else float(self.rt[hi])
        return max(0.0, right - left)

    @property
    def fwhm_is_truncated(self) -> bool:
        """True when the half-maximum region reaches an end of the trace.

        The peak is cut off by the extraction window, so `fwhm` and
        `points_across_peak` are lower bounds, not measurements.
        """
        span = self._half_max_span()
        if span is None:
            return False
        _apex, _half, lo, hi = span
        return lo == 0 or hi == len(self.intensity) - 1

    @property
    def points_across_peak(self) -> int:
        """Scans acquired within the full width at half maximum.

        This is the quantity instrument methods expose as "desired minimum points
        across the peak"; comparing the measured value against that setting is how
        you tell whether the chromatography is being sampled adequately, and whether
        a longer gradient has room to broaden peaks further.
        """
        span = self._half_max_span()
        if span is None:
            return 0
        _apex, _half, lo, hi = span
        return int(hi - lo + 1)


@dataclass
class ChromatogramSet:
    """A shared retention-time axis plus one or more co-extracted XICs.

    Produced by :func:`extract_xics` in a single pass so every trace (and the
    TIC / base-peak lanes) share the same scan grid -- ready to feed a stacked
    multi-lane chromatogram panel.
    """

    rt: np.ndarray
    tic: np.ndarray
    base_peak: np.ndarray
    xics: Dict[str, XIC]
    ms_level: int = 0
    scan_nums: Optional[np.ndarray] = None

    def __getitem__(self, key: str) -> XIC:
        return self.xics[key]

    def __iter__(self):
        return iter(self.xics.values())

    def __len__(self) -> int:
        return len(self.xics)

    def names(self) -> List[str]:
        return list(self.xics.keys())

    @property
    def n_scans(self) -> int:
        return int(len(self.rt))


def _normalize_targets(targets: TargetSpec) -> Tuple[List[str], List[float]]:
    """Accept a {name: mz} mapping, a list of (name, mz) pairs, or a list of
    bare m/z floats (auto-named by value)."""
    if isinstance(targets, Mapping):
        return list(targets.keys()), [float(v) for v in targets.values()]
    names: List[str] = []
    mzs: List[float] = []
    for item in targets:
        if isinstance(item, (tuple, list)) and len(item) == 2:
            names.append(str(item[0]))
            mzs.append(float(item[1]))
        else:
            mz = float(item)  # type: ignore[arg-type]
            mzs.append(mz)
            names.append(f"{mz:.4f}")
    return names, mzs


def _window_intensity(mz: np.ndarray, intensity: np.ndarray, target_mz: float,
                      tolerance: float, unit: str, mode: str) -> float:
    """Intensity of *target_mz* within tolerance for one spectrum.

    ``mode='sum'`` sums every peak inside the window (standard XIC behaviour,
    robust to profile/split peaks); ``mode='max'`` takes the apex peak.
    """
    if len(mz) == 0:
        return 0.0
    if unit == "ppm":
        tol_da = target_mz * tolerance / 1e6
    elif unit == "Da":
        tol_da = tolerance
    else:
        raise ValueError(f"Unknown unit '{unit}'. Use 'ppm' or 'Da'.")
    mask = np.abs(mz - target_mz) <= tol_da
    if not mask.any():
        return 0.0
    vals = intensity[mask]
    if mode == "sum":
        return float(vals.sum())
    elif mode == "max":
        return float(vals.max())
    raise ValueError(f"Unknown mode '{mode}'. Use 'sum' or 'max'.")


def _activation_match(spec, activation: str) -> bool:
    act = activation.lower()
    at = getattr(spec, "activation_type", "") or ""
    if at.lower() == act:
        return True
    fs = (getattr(spec, "filter_string", "") or "").lower()
    return act in fs


def extract_xics(source, targets: TargetSpec, *,
                 tolerance: float = 20.0,
                 unit: str = "ppm",
                 ms_level: Optional[int] = 1,
                 rt_range: Optional[Tuple[float, float]] = None,
                 precursor_mz: Optional[float] = None,
                 precursor_tol: float = 0.7,
                 mode: str = "sum",
                 activation: Optional[str] = None,
                 collect_tic: bool = True) -> ChromatogramSet:
    """Extract several XICs plus the TIC / base-peak traces in one pass.

    Args:
        source: An open reader (anything exposing ``iter_spectra()`` -- an
            ``MzMLReader`` or ``SpectrumCache``) or a path/str. A path is
            opened cache-aware via ``open_spectra`` and closed on exit.
        targets: Target ions -- a ``{name: mz}`` mapping, ``(name, mz)`` pairs,
            or bare m/z floats.
        tolerance: Mass tolerance for the extraction window.
        unit: ``'ppm'`` (default) or ``'Da'``.
        ms_level: Only scans of this MS level contribute (``1`` for precursor
            XICs, ``2`` for fragment/oxonium XICs). ``None`` uses every scan.
        rt_range: ``(min, max)`` retention-time window (reader's native unit,
            usually minutes). ``None`` uses the whole run.
        precursor_mz: If given, keep only scans whose precursor m/z is within
            ``precursor_tol`` Da (useful for a single-precursor fragment XIC).
        precursor_tol: Precursor match tolerance in Da (default 0.7).
        mode: ``'sum'`` (default) or ``'max'`` intensity within the window.
        activation: If given (e.g. ``'HCD'``), keep only scans of that
            activation type -- skips ETD/EThcD scans in hybrid runs.
        collect_tic: Also collect per-scan TIC and base-peak traces.

    Returns:
        A :class:`ChromatogramSet` with a shared RT axis, the TIC / base-peak
        traces, and one :class:`XIC` per target (keyed by name).
    """
    names, mzs = _normalize_targets(targets)

    rts: List[float] = []
    scans: List[int] = []
    tics: List[float] = []
    bps: List[float] = []
    lanes: List[List[float]] = [[] for _ in mzs]

    def _run(reader):
        for spec in reader.iter_spectra():
            if ms_level is not None and spec.ms_level != ms_level:
                continue
            if activation is not None and not _activation_match(spec, activation):
                continue
            if rt_range is not None and not (rt_range[0] <= spec.rt <= rt_range[1]):
                continue
            if precursor_mz is not None:
                pmz = getattr(spec, "precursor_mz", 0.0) or 0.0
                if pmz <= 0 or abs(pmz - precursor_mz) > precursor_tol:
                    continue
            rts.append(spec.rt)
            scans.append(spec.scan_num)
            mz_a = spec.mz
            int_a = spec.intensity
            if collect_tic:
                tic = spec.tic if spec.tic else (float(int_a.sum()) if len(int_a) else 0.0)
                bp = (spec.base_peak_intensity if spec.base_peak_intensity
                      else (float(int_a.max()) if len(int_a) else 0.0))
                tics.append(tic)
                bps.append(bp)
            for i, t in enumerate(mzs):
                lanes[i].append(_window_intensity(mz_a, int_a, t, tolerance, unit, mode))

    # Reader passed in -> use directly (do not close the caller's reader).
    # Path passed in -> open cache-aware and close on exit.
    if hasattr(source, "iter_spectra"):
        _run(source)
    else:
        from .spectrum_cache import open_spectra
        with open_spectra(str(source)) as reader:
            _run(reader)

    rt = np.asarray(rts, dtype=float)
    order = np.argsort(rt, kind="stable") if rt.size else np.array([], dtype=int)
    rt = rt[order]
    scan_arr = np.asarray(scans)[order] if scans else np.array([], dtype=int)
    tic_arr = np.asarray(tics, dtype=float)[order] if collect_tic and tics else np.array([], dtype=float)
    bp_arr = np.asarray(bps, dtype=float)[order] if collect_tic and bps else np.array([], dtype=float)

    xics: Dict[str, XIC] = {}
    for i, (nm, t) in enumerate(zip(names, mzs)):
        inten = np.asarray(lanes[i], dtype=float)[order] if lanes[i] else np.array([], dtype=float)
        xics[nm] = XIC(target_mz=t, rt=rt, intensity=inten, name=nm,
                       ms_level=ms_level if ms_level is not None else 0)

    return ChromatogramSet(rt=rt, tic=tic_arr, base_peak=bp_arr, xics=xics,
                           ms_level=ms_level if ms_level is not None else 0,
                           scan_nums=scan_arr)


def extract_xic(source, target_mz: float, *,
                name: str = "",
                tolerance: float = 20.0,
                unit: str = "ppm",
                ms_level: Optional[int] = 1,
                rt_range: Optional[Tuple[float, float]] = None,
                precursor_mz: Optional[float] = None,
                precursor_tol: float = 0.7,
                mode: str = "sum",
                activation: Optional[str] = None) -> XIC:
    """Extract a single XIC for one target m/z. See :func:`extract_xics` for
    argument semantics. Returns one :class:`XIC`."""
    nm = name or f"{float(target_mz):.4f}"
    cs = extract_xics(source, {nm: target_mz}, tolerance=tolerance, unit=unit,
                      ms_level=ms_level, rt_range=rt_range, precursor_mz=precursor_mz,
                      precursor_tol=precursor_tol, mode=mode, activation=activation,
                      collect_tic=False)
    return cs.xics[nm]
