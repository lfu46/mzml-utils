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

import inspect
import math
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from .constants import NEUTRON_MASS, PROTON

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
        """Trapezoidal area over the WHOLE trace (every point in the extraction window).

        This is not a peak area: a second elution or background inside the window is
        integrated too. For the area between the peak's boundaries use
        :func:`peak_metrics`.
        """
        if len(self.rt) < 2:
            return 0.0
        return float(_trapz(self.intensity, self.rt))

    def _half_max_span(self, apex: Optional[int] = None):
        """(apex_idx, half, lo_idx, hi_idx) of the CONTIGUOUS half-maximum region.

        Contiguous on purpose: an XIC often contains more than one peak (a second
        elution, or an interfering species at the same m/z), and counting every point
        above half maximum anywhere in the trace would merge them into one absurdly
        wide peak. Only the region containing the apex is the peak.

        ``apex`` is the index of the peak to measure; ``None`` takes the tallest point.
        """
        n = len(self.intensity)
        if n == 0:
            return None
        if apex is None:
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

    def _half_max_crossings(self, apex: Optional[int] = None):
        """(left_rt, right_rt, lo_idx, hi_idx) where the trace crosses half maximum.

        The crossings are linearly interpolated between the bracketing scans. At an
        end of the trace the crossing is snapped to that end point (a lower bound).
        """
        span = self._half_max_span(apex)
        if span is None or len(self.rt) < 2:
            return None
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
        return left, right, lo, hi

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
        cr = self._half_max_crossings()
        if cr is None:
            return 0.0
        left, right, _lo, _hi = cr
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
        # Push the MS-level / RT filters into the reader when it accepts them (both
        # built-in readers do), so a SpectrumCache never decodes the scans it would
        # drop. The Python-side checks below stay as the guard for any other reader.
        accepted = inspect.signature(reader.iter_spectra).parameters
        pushed = {k: v for k, v in (("ms_level", ms_level), ("rt_range", rt_range))
                  if v is not None and k in accepted}
        for spec in reader.iter_spectra(**pushed):
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


# ---------------------------------------------------------------------------
# Peak boundaries and peak-quality metrics
# ---------------------------------------------------------------------------
#
# Ported from the R sources, not their documentation, so the numbers match:
#   * boundaries  -- Chromatograms (Bioconductor, Louail et al. Anal Chem 2026),
#                    R/helpers.R `.peak_boundary_one`
#   * snr, prominence -- MsQuality R/function_Chromatograms_metrics.R
#                    (`signalToNoiseRatio`, `.peakProminenceSingle`)
#   * beta_cor, beta_snr -- MetaboCoreUtils R/peak-shape-quality.R `betaValues`
#                    (Kumler, Hazelton & Ingalls, BMC Bioinformatics 2023, 24:404)
# tests/test_xic.py pins the port against values produced by running those R
# functions verbatim.

_MAD_CONSTANT = 1.4826  # R stats::mad default: consistent with the SD of a normal


def _bounds_from_apex(inten: np.ndarray, apex: int, baseline: float, threshold: float,
                      baseline_threshold: float, nearest_valley: bool) -> Tuple[int, int]:
    """Boundary indices of the peak at ``apex`` (steps 2-4 of :func:`peak_boundary`)."""
    n = len(inten)
    height = float(inten[apex]) - baseline
    valley_ok = baseline + height * baseline_threshold
    if nearest_valley:
        # Leave a flat top, then descend while strictly falling: the first point whose
        # outward neighbour is not lower is the valley's nearest point.
        left = apex
        while left > 0 and inten[left - 1] == inten[apex]:
            left -= 1
        while left > 0 and inten[left - 1] < inten[left]:
            left -= 1
        right = apex
        while right < n - 1 and inten[right + 1] == inten[apex]:
            right += 1
        while right < n - 1 and inten[right + 1] < inten[right]:
            right += 1
    else:
        left = apex
        while left > 0 and not (inten[left] < inten[left - 1]
                                and (left == n - 1 or inten[left] <= inten[left + 1])):
            left -= 1
        right = apex
        while right < n - 1 and not ((right == 0 or inten[right] < inten[right - 1])
                                     and inten[right] <= inten[right + 1]):
            right += 1

    if not (inten[left] <= valley_ok and inten[right] <= valley_ok):
        thresh = baseline + height * threshold
        left_cand = np.flatnonzero(inten[:apex + 1] <= thresh)
        right_cand = np.flatnonzero(inten[apex:] <= thresh)
        left = int(left_cand[-1]) if left_cand.size else 0
        right = apex + int(right_cand[0]) if right_cand.size else n - 1
    return int(left), int(right)


def peak_boundary(xic: "XIC", *,
                  contains_rt: Optional[float] = None,
                  apex_window: Optional[Tuple[float, float]] = None,
                  threshold: float = 0.1,
                  baseline_threshold: float = 0.1,
                  baseline_quantile: float = 0.1,
                  nearest_valley: bool = False) -> Optional[Tuple[int, int, int]]:
    """Left boundary, apex and right boundary of one chromatographic peak, as indices.

    The algorithm of ``Chromatograms::peakBoundary``:

    1. baseline = the ``baseline_quantile`` quantile of the trace's intensities;
    2. walk outwards from the apex to the first local minimum (valley) on each side;
    3. keep the two valleys if BOTH are at or below
       ``baseline + baseline_threshold * (apex - baseline)``;
    4. otherwise -- for either side failing, both sides are redone -- take the nearest
       points on each side at or below ``baseline + threshold * (apex - baseline)``,
       or the trace ends.

    With the defaults this reproduces R exactly: the apex is the tallest point in the
    trace. Three options deviate from R, all off by default.

    ``contains_rt`` selects the peak that CONTAINS a retention time, e.g. the PSM's MS2
    scan time. A glycopeptide XIC often holds more than one peak (an isomer, a +/-1 Da
    neighbour, a second elution), and the tallest need not be the assigned species.
    Peaks are taken tallest first; each is bounded as in steps 2-4 on the whole trace,
    and the first whose boundaries enclose ``contains_rt`` is returned, with its apex
    re-taken as the tallest point inside those boundaries. A peak that does not enclose
    it is set aside and the next tallest remaining point is tried. When ``contains_rt``
    lies on a tail below ``threshold`` of the main apex, R's boundaries exclude it;
    the enclosing boundaries then come from a tail point's own, lower threshold and are
    wider than R would draw for the apex they are reported with. Selecting by
    containment rather than by a window around the MS2 time matters in DDA: the MS2
    can be triggered on a peak's tail, a minute or more from its apex, and a window
    would then measure a point on the tail as if it were the apex.

    ``apex_window=(lo, hi)`` takes the apex as the tallest point inside the window.

    ``nearest_valley=True`` changes the walk on a FLAT valley. R's walk is asymmetric
    there: the right side stops at the first flat point, but the left side walks
    through the whole flat run (it keeps a flat minimum at its first, leftmost index).
    A peak on a perfectly flat floor -- in particular a zero-filled XIC, where every
    scan without the ion is exactly 0 -- therefore gets its left boundary at the START
    OF THE TRACE, and every in-boundary metric is computed over the empty floor (a
    perfect Gaussian on a flat floor scores beta_cor 0.35 in R). With
    ``nearest_valley=True`` both sides stop at the flat run's point nearest the apex.
    Traces with a noisy baseline rarely have exact ties and give the same result
    either way.

    Any valley walk, R's included, stops at an internal DROPOUT: a scan inside the
    peak where the ion is absent (intensity 0) is a valley at the baseline. Check the
    trace before reading a boundary as the end of elution.

    Returns ``None`` when there are fewer than 3 points, the apex intensity is zero,
    the apex window holds no point, or no peak above baseline contains
    ``contains_rt``. The trace must not contain NaN; extraction writes 0 for an absent
    ion.
    """
    if contains_rt is not None and apex_window is not None:
        raise ValueError("peak_boundary: give contains_rt or apex_window, not both")
    inten = np.asarray(xic.intensity, dtype=float)
    rt = np.asarray(xic.rt, dtype=float)
    n = len(inten)
    if np.isnan(inten).any():
        raise ValueError("peak_boundary: intensity contains NaN (write 0 for an absent ion)")
    if n < 3:
        return None
    baseline = float(np.quantile(inten, baseline_quantile))  # numpy 'linear' == R type 7
    args = (baseline, threshold, baseline_threshold, nearest_valley)

    if contains_rt is not None:
        t = float(contains_rt)
        if not (rt[0] <= t <= rt[-1]):
            return None
        free = np.ones(n, dtype=bool)
        while True:
            cand = np.flatnonzero(free & (inten > baseline))
            if cand.size == 0:
                return None
            apex = int(cand[int(np.argmax(inten[cand]))])
            left, right = _bounds_from_apex(inten, apex, *args)
            if rt[left] <= t <= rt[right]:
                # The apex is the tallest point inside the accepted boundaries. It differs
                # from `apex` when `t` lies on a tail below the main peak's threshold: the
                # tail point's own (lower) threshold then draws boundaries that span the
                # main peak as well.
                top = left + int(np.argmax(inten[left:right + 1]))
                return left, top, right
            free[left:right + 1] = False

    if apex_window is None:
        apex = int(np.argmax(inten))
    else:
        inside = np.flatnonzero((rt >= apex_window[0]) & (rt <= apex_window[1]))
        if inside.size == 0:
            return None
        apex = int(inside[int(np.argmax(inten[inside]))])
        if inten[apex] <= baseline:
            return None
    if inten[apex] == 0.0:
        return None
    left, right = _bounds_from_apex(inten, apex, *args)
    return left, apex, right


def _scale_zero_one(v: np.ndarray) -> np.ndarray:
    return (v - v.min()) / (v.max() - v.min())


def _beta_pdf(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Beta(a, b) density on [0, 1] (a, b > 1, so it is 0 at both ends)."""
    log_norm = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    out = np.zeros_like(x, dtype=float)
    inner = (x > 0.0) & (x < 1.0)
    xi = x[inner]
    out[inner] = np.exp(log_norm + (a - 1.0) * np.log(xi) + (b - 1.0) * np.log1p(-xi))
    return out


def beta_values(intensity, rtime=None,
                skews: Sequence[float] = (3.0, 3.5, 4.0, 4.5, 5.0)) -> Tuple[float, float]:
    """Kumler's beta peak-shape metrics: ``(beta_cor, beta_snr)``.

    A port of ``MetaboCoreUtils::betaValues``. Zero intensities are dropped, retention
    time is rescaled to [0, 1], and the trace is correlated (Pearson) with Beta(skew, 5)
    densities for each ``skew`` (values below 5 are right-skewed, i.e. tailing).

    * ``beta_cor`` -- the best correlation. Near 1 for a clean single peak; low for
      noise or a trace with several maxima.
    * ``beta_snr`` -- ``log10(max(intensity) / sd(diff(scaled best curve - scaled
      intensity)))``. HIGHER is better. MsQuality labels this column
      ``gaussian_residuals`` and documents it as a residual SD where lower is better;
      its code returns this log10 ratio. Note the numerator is the RAW maximum
      intensity while the noise is computed on 0-1 scaled curves, so the value grows
      with absolute signal and is not a pure shape score.

    Returns ``(nan, nan)`` with fewer than 5 non-zero points (R returns NA), or when
    the non-zero intensities are all equal (R's correlation is NA there).
    """
    y = np.asarray(intensity, dtype=float)
    x = np.arange(1, len(y) + 1, dtype=float) if rtime is None else np.asarray(rtime, dtype=float)
    keep = y > 0
    y, x = y[keep], x[keep]
    if len(y) < 5 or y.max() == y.min():
        return float("nan"), float("nan")
    xs = _scale_zero_one(x)
    curves = [_beta_pdf(xs, float(s), 5.0) for s in skews]
    cors = [float(np.corrcoef(y, c)[0, 1]) for c in curves]
    best = int(np.argmax(cors))  # first maximum, as R's which.max
    noise = float(np.std(np.diff(_scale_zero_one(curves[best]) - _scale_zero_one(y)), ddof=1))
    with np.errstate(divide="ignore"):
        beta_snr = float(np.log10(y.max() / noise)) if noise > 0 else float("inf")
    return cors[best], beta_snr


@dataclass
class PeakMetrics:
    """Boundaries and quality metrics of one chromatographic peak (see :func:`peak_metrics`).

    Retention times are in the trace's unit (minutes). NaN means "not determinable";
    a metric that needs a peak is NaN when no boundary was found.
    """

    apex_rt: float = float("nan")
    apex_intensity: float = 0.0
    left_rt: float = float("nan")
    right_rt: float = float("nan")
    boundary_is_truncated: bool = False
    """A boundary sits on an end of the trace: the extraction window cut the peak."""
    n_points: int = 0
    """Scans between the boundaries, inclusive."""
    area: float = 0.0
    """Trapezoidal area between the boundaries."""
    fwhm: float = float("nan")
    fwhm_is_truncated: bool = False
    points_across_peak: int = 0
    asymmetry: float = float("nan")
    """(right half-max crossing - apex) / (apex - left half-max crossing). 1 is
    symmetric; above 1 tails, below 1 fronts."""
    snr: float = float("nan")
    """Apex intensity / MAD of the whole trace (MsQuality ``signalToNoiseRatio``)."""
    prominence: float = float("nan")
    """(apex - q10) / q10 over the whole trace (MsQuality ``peakProminence``)."""
    beta_cor: float = float("nan")
    beta_snr: float = float("nan")

    @property
    def found(self) -> bool:
        return self.apex_intensity > 0.0 and not math.isnan(self.left_rt)

    @property
    def width(self) -> float:
        return self.right_rt - self.left_rt

    def as_dict(self) -> Dict[str, object]:
        d = {f: getattr(self, f) for f in self.__dataclass_fields__}
        d["width"] = self.width
        return d


def peak_metrics(xic: "XIC", *,
                 contains_rt: Optional[float] = None,
                 apex_window: Optional[Tuple[float, float]] = None,
                 threshold: float = 0.1,
                 baseline_threshold: float = 0.1,
                 baseline_quantile: float = 0.1,
                 nearest_valley: bool = False) -> PeakMetrics:
    """Boundaries, in-boundary area and quality metrics of one peak in an XIC.

    The boundary arguments are those of :func:`peak_boundary`: ``contains_rt`` (or
    ``apex_window``) selects which peak is measured, and ``nearest_valley=True`` is
    recommended for zero-filled traces (see there). Metric definitions follow their R sources:

    * ``area``, ``n_points`` and the ``beta_*`` metrics use the points between the
      boundaries (MsQuality ``gaussianSimilarity`` does the same);
    * ``fwhm``, ``points_across_peak``, ``asymmetry`` use the contiguous half-maximum
      region around the apex, crossings interpolated (as :attr:`XIC.fwhm`);
    * ``snr`` and ``prominence`` use the whole trace, as MsQuality does by default,
      with the APEX intensity as numerator. With neither selector the apex is the
      trace maximum, so both equal MsQuality's values. Both are NaN when more than
      half (MAD) or a tenth (q10) of the trace is zero, as MsQuality returns NA --
      common for a zero-filled XIC over a wide window.

    These are measurements, not verdicts: no threshold is applied here.
    """
    b = peak_boundary(xic, contains_rt=contains_rt, apex_window=apex_window, threshold=threshold,
                      baseline_threshold=baseline_threshold,
                      baseline_quantile=baseline_quantile,
                      nearest_valley=nearest_valley)
    if b is None:
        return PeakMetrics()
    left, apex, right = b
    rt = np.asarray(xic.rt, dtype=float)
    inten = np.asarray(xic.intensity, dtype=float)
    apex_int = float(inten[apex])
    apex_rt = float(rt[apex])

    seg_rt, seg_int = rt[left:right + 1], inten[left:right + 1]
    area = float(_trapz(seg_int, seg_rt)) if right > left else 0.0

    fwhm = asym = float("nan")
    fwhm_trunc, pts = False, 0
    cr = xic._half_max_crossings(apex)
    if cr is not None:
        l_x, r_x, lo, hi = cr
        fwhm = max(0.0, r_x - l_x)
        fwhm_trunc = lo == 0 or hi == len(inten) - 1
        pts = int(hi - lo + 1)
        if apex_rt - l_x > 0 and r_x - apex_rt > 0:
            asym = (r_x - apex_rt) / (apex_rt - l_x)

    med = float(np.median(inten))
    mad = _MAD_CONSTANT * float(np.median(np.abs(inten - med)))
    snr = apex_int / mad if mad > 0 else float("nan")
    q10 = float(np.quantile(inten, baseline_quantile))
    prominence = (apex_int - q10) / q10 if q10 > 0 else float("nan")

    beta_cor, beta_snr = beta_values(seg_int, seg_rt) if len(seg_int) >= 5 else (float("nan"),) * 2

    return PeakMetrics(apex_rt=apex_rt, apex_intensity=apex_int,
                       left_rt=float(rt[left]), right_rt=float(rt[right]),
                       boundary_is_truncated=(left == 0 or right == len(inten) - 1),
                       n_points=int(right - left + 1), area=area,
                       fwhm=fwhm, fwhm_is_truncated=fwhm_trunc, points_across_peak=pts,
                       asymmetry=asym, snr=snr, prominence=prominence,
                       beta_cor=beta_cor, beta_snr=beta_snr)


# ---------------------------------------------------------------------------
# Isotope co-elution -- trace-level precursor evidence
# ---------------------------------------------------------------------------

def coelution_score(a: XIC, b: XIC,
                    rt_range: Optional[Tuple[float, float]] = None) -> float:
    """Cosine similarity of two XICs that share one retention-time axis.

    1.0 means the traces rise and fall together (isotopes of one species); near 0
    means they do not co-elute. The dot product and both norms are taken over the
    SAME points, so a trace that only partly overlaps the other is scored for what
    it is rather than penalised twice.

    The traces must come from one :func:`extract_xics` call (same scan grid).
    Returns 0.0 when either trace is empty or flat zero in the window.
    """
    if len(a.intensity) != len(b.intensity):
        raise ValueError("XICs must share one RT axis (extract them in one extract_xics call)")
    x, y = a.intensity, b.intensity
    if rt_range is not None:
        keep = (a.rt >= rt_range[0]) & (a.rt <= rt_range[1])
        x, y = x[keep], y[keep]
    norm = float(np.linalg.norm(x) * np.linalg.norm(y))
    return float(np.dot(x, y) / norm) if norm > 0 else 0.0


@dataclass
class IsotopeTraceEvidence:
    """MS1 trace-level evidence for one precursor. Keys are isotope indices
    relative to the assigned monoisotopic peak: -1, 0, 1, 2, ..."""

    mono_mz: float
    charge: int
    apex_rt: float
    """Apex of the monoisotopic trace inside the window (0.0 when it is empty)."""
    coelution: Dict[int, float]
    """Cosine of each isotope trace with the monoisotopic trace (index 0 omitted)."""
    peak_intensity: Dict[int, float]
    """Intensity of each isotope summed over the monoisotopic peak's half-maximum span."""
    pattern_cosine: float
    """Cosine of the observed M, M+1, ... intensities against the theoretical envelope."""
    m_minus_1_coelutes: bool
    """True when a trace one isotope step BELOW the assigned monoisotopic m/z co-elutes
    with it -- the assigned peak is then probably M+1 of a lighter species."""
    chromatograms: ChromatogramSet

    @property
    def found(self) -> bool:
        return self.peak_intensity.get(0, 0.0) > 0.0


def isotope_trace_evidence(source, mono_mz: float, charge: int, *,
                           rt_center: float, rt_halfwidth: float = 1.0,
                           n_isotopes: int = 3,
                           tolerance: float = 10.0, unit: str = "ppm",
                           distribution=None,
                           min_coelution: float = 0.7,
                           min_m_minus_1_ratio: float = 0.1) -> IsotopeTraceEvidence:
    """Do the isotopes of a precursor co-elute, and is the monoisotopic peak the right one?

    A single MS1 scan cannot tell an isotope from an unrelated ion that happens to
    sit one neutron away; a chromatographic trace can, because isotopes of one
    species share an elution profile. One :func:`extract_xics` pass pulls the traces
    for M-1, M, M+1, ... inside ``rt_center +/- rt_halfwidth`` and scores them.

    Isotope m/z values are anchored on ``mono_mz`` (pass the OBSERVED value when you
    have it, so calibration error cancels in the spacing). The spacing is
    ``NEUTRON_MASS / charge``, or the centroid spacing of ``distribution`` when one is
    given -- needed for a composition whose envelope is not averagine-like (a halogen
    tag, a heavy label). The M-1 trace is always one ``NEUTRON_MASS / charge`` below.

    Args:
        source: An open reader or a path (see :func:`extract_xics`).
        mono_mz: Assigned monoisotopic m/z.
        charge: Precursor charge.
        rt_center, rt_halfwidth: Retention window, in the reader's unit (minutes).
        n_isotopes: Peaks of the envelope to trace (M through M+(n-1)).
        tolerance, unit: Extraction window per trace.
        distribution: Optional :class:`~mzml_utils.isotopes.IsotopeDistribution` of the
            neutral composition. ``None`` uses the averagine of ``deisotope``.
        min_coelution: Cosine an M-1 trace must reach to count as co-eluting.
        min_m_minus_1_ratio: Minimum M-1 / M intensity ratio for the same verdict; a
            real lighter monoisotopic peak is rarely below this under ~8 kDa. Both
            are screening defaults, not calibrated thresholds -- read the numbers.
    """
    if charge < 1:
        raise ValueError("charge must be >= 1")
    step = NEUTRON_MASS / charge
    if distribution is not None:
        n_isotopes = min(n_isotopes, distribution.n_peaks)
        offsets = [float(distribution.spacing[k]) / charge for k in range(n_isotopes)]
        theoretical = np.asarray(distribution.abundance[:n_isotopes], dtype=float)
    else:
        from .deisotope import _poisson_distribution
        offsets = [k * step for k in range(n_isotopes)]
        theoretical = _poisson_distribution((mono_mz - PROTON) * charge, n_isotopes)

    targets = {"M-1": mono_mz - step}
    targets.update({f"M+{k}": mono_mz + offsets[k] for k in range(n_isotopes)})
    cs = extract_xics(source, targets, tolerance=tolerance, unit=unit, ms_level=1,
                      rt_range=(rt_center - rt_halfwidth, rt_center + rt_halfwidth),
                      collect_tic=False)

    mono = cs["M+0"]
    index = {-1: "M-1", **{k: f"M+{k}" for k in range(n_isotopes)}}
    span = mono._half_max_span()
    if span is None:
        lo, hi = 0, -1
    else:
        _apex, _half, lo, hi = span
    peak_intensity = {k: float(cs[nm].intensity[lo:hi + 1].sum()) for k, nm in index.items()}
    coelution = {k: coelution_score(mono, cs[nm]) for k, nm in index.items() if k != 0}

    observed = np.array([peak_intensity[k] for k in range(n_isotopes)])
    norm = float(np.linalg.norm(observed) * np.linalg.norm(theoretical))
    pattern_cosine = float(np.dot(observed, theoretical) / norm) if norm > 0 else 0.0

    m0 = peak_intensity[0]
    m_minus_1 = (m0 > 0 and coelution[-1] >= min_coelution
                 and peak_intensity[-1] / m0 >= min_m_minus_1_ratio)

    return IsotopeTraceEvidence(
        mono_mz=float(mono_mz), charge=int(charge),
        apex_rt=mono.apex_rt if span is not None else 0.0,
        coelution=coelution, peak_intensity=peak_intensity,
        pattern_cosine=pattern_cosine, m_minus_1_coelutes=bool(m_minus_1),
        chromatograms=cs)
