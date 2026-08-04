"""
Spectral averaging with outlier rejection, weighting, and normalization.

Combines multiple spectra into a single averaged spectrum using configurable
binning, outlier rejection, weighting, and normalization strategies.

Reference: mzLib SpectralAveraging (Smith lab, UW-Madison).
Outlier rejection methods adapted from PixInsight ImageIntegration.
"""

import numpy as np
from enum import Enum
from typing import List, Tuple, Optional


class OutlierRejectionMethod(Enum):
    NONE = "none"
    MIN_MAX = "min_max"
    PERCENTILE = "percentile"
    SIGMA_CLIPPING = "sigma_clipping"
    WINSORIZED_SIGMA = "winsorized_sigma"
    AVERAGED_SIGMA = "averaged_sigma"
    BELOW_THRESHOLD = "below_threshold"


class WeightingScheme(Enum):
    EVEN = "even"
    TIC = "tic"


class AvgNormalization(Enum):
    NONE = "none"
    ABSOLUTE_TIC = "absolute_tic"
    RELATIVE_TIC = "relative_tic"
    RELATIVE_INTENSITY = "relative_intensity"


# --- Outlier rejection implementations ---

def _reject_none(values: np.ndarray, **kwargs) -> np.ndarray:
    return values


def _reject_min_max(values: np.ndarray, **kwargs) -> np.ndarray:
    if len(values) <= 2:
        return values
    mn, mx = values.min(), values.max()
    return values[(values > mn) & (values < mx)]


def _reject_percentile(values: np.ndarray, percentile: float = 0.9, **kwargs) -> np.ndarray:
    if len(values) <= 2:
        return values
    trim = (1.0 - percentile) / 2.0
    med = np.median(values)
    high_pct = 1.0 - trim
    high_cut = med * (1.0 + high_pct)
    low_cut = med * (1.0 - high_pct)
    return values[(values > low_cut) & (values < high_cut)]


def _sigma_clip_once(values: np.ndarray, sigma_low: float, sigma_high: float) -> np.ndarray:
    med = np.median(values)
    std = np.std(values)
    if std == 0:
        return values
    keep = np.ones(len(values), dtype=bool)
    keep[(med - values) / std > sigma_low] = False
    keep[(values - med) / std > sigma_high] = False
    return values[keep]


def _reject_sigma_clipping(values: np.ndarray, sigma_low: float = 3.0,
                            sigma_high: float = 3.0, **kwargs) -> np.ndarray:
    prev_len = len(values) + 1
    result = values.copy()
    while len(result) < prev_len and len(result) > 1:
        prev_len = len(result)
        result = _sigma_clip_once(result, sigma_low, sigma_high)
    return result


def _reject_winsorized_sigma(values: np.ndarray, sigma_low: float = 3.0,
                              sigma_high: float = 3.0, **kwargs) -> np.ndarray:
    """Winsorized sigma clipping with Huber loop for robust sigma estimation."""
    if len(values) <= 2:
        return values
    result = values.copy()
    prev_len = len(result) + 1

    while len(result) < prev_len and len(result) > 1:
        prev_len = len(result)
        nonzero = result[result > 0]
        if len(nonzero) == 0:
            break

        med = np.median(nonzero)
        std = np.std(nonzero)
        if std == 0:
            break

        # Huber loop: converge on robust sigma
        for _ in range(100):
            winsorized = np.clip(result, med - 1.5 * std, med + 1.5 * std)
            med_new = np.median(winsorized)
            std_new = np.std(winsorized) * 1.134  # PixInsight correction
            if std_new == 0 or abs(std_new - std) / std < 0.00005:
                std = std_new if std_new > 0 else std
                med = med_new
                break
            std = std_new
            med = med_new

        if std == 0:
            break
        result = _sigma_clip_once(result, sigma_low, sigma_high)

    return result


def _reject_averaged_sigma(values: np.ndarray, sigma_low: float = 3.0,
                            sigma_high: float = 3.0, **kwargs) -> np.ndarray:
    """Averaged sigma clipping where sigma depends on sqrt(median)."""
    if len(values) <= 2:
        return values
    result = values.copy()
    prev_len = len(result) + 1
    med = np.median(result)
    if med <= 0:
        return result
    s = np.sqrt(np.sum((result - med) ** 2 / med) / max(len(result) - 1, 1))

    while len(result) < prev_len and len(result) > 1:
        prev_len = len(result)
        nonzero = result[result > 0]
        if len(nonzero) == 0:
            break
        med = np.median(nonzero)
        if med <= 0:
            break
        sigma = s * np.sqrt(med)
        if sigma == 0:
            break
        keep = np.ones(len(result), dtype=bool)
        keep[(med - result) / sigma > sigma_low] = False
        keep[(result - med) / sigma > sigma_high] = False
        result = result[keep]

    return result


def _reject_below_threshold(values: np.ndarray, cutoff: float = 0.2, **kwargs) -> np.ndarray:
    """If fraction of non-zero values <= cutoff, zero out entire bin."""
    nonzero_frac = np.count_nonzero(values) / len(values) if len(values) > 0 else 0
    if nonzero_frac <= cutoff:
        return np.zeros_like(values)
    return values


_REJECTION_FUNCS = {
    OutlierRejectionMethod.NONE: _reject_none,
    OutlierRejectionMethod.MIN_MAX: _reject_min_max,
    OutlierRejectionMethod.PERCENTILE: _reject_percentile,
    OutlierRejectionMethod.SIGMA_CLIPPING: _reject_sigma_clipping,
    OutlierRejectionMethod.WINSORIZED_SIGMA: _reject_winsorized_sigma,
    OutlierRejectionMethod.AVERAGED_SIGMA: _reject_averaged_sigma,
    OutlierRejectionMethod.BELOW_THRESHOLD: _reject_below_threshold,
}


# --- Normalization ---

def _normalize_spectra(intensity_arrays: List[np.ndarray],
                       method: AvgNormalization) -> Tuple[List[np.ndarray], float]:
    """Normalize spectra. Returns (normalized_arrays, average_tic)."""
    tics = [y.sum() for y in intensity_arrays]
    avg_tic = np.mean(tics) if tics else 1.0

    if method == AvgNormalization.NONE:
        return [y.copy() for y in intensity_arrays], avg_tic

    elif method == AvgNormalization.ABSOLUTE_TIC:
        return [y / t if t > 0 else y.copy()
                for y, t in zip(intensity_arrays, tics)], avg_tic

    elif method == AvgNormalization.RELATIVE_TIC:
        return [y / t * avg_tic if t > 0 else y.copy()
                for y, t in zip(intensity_arrays, tics)], avg_tic

    elif method == AvgNormalization.RELATIVE_INTENSITY:
        return [y / y.max() if y.max() > 0 else y.copy()
                for y in intensity_arrays], avg_tic

    return [y.copy() for y in intensity_arrays], avg_tic


# --- Weighting ---

def _compute_weights(intensity_arrays: List[np.ndarray],
                     scheme: WeightingScheme) -> np.ndarray:
    """Compute per-spectrum weights."""
    n = len(intensity_arrays)
    if scheme == WeightingScheme.EVEN:
        return np.ones(n)
    elif scheme == WeightingScheme.TIC:
        tics = np.array([y.sum() for y in intensity_arrays])
        mx = tics.max()
        return tics / mx if mx > 0 else np.ones(n)
    return np.ones(n)


# --- Main averaging function ---

def average_spectra(mz_arrays: List[np.ndarray],
                    intensity_arrays: List[np.ndarray],
                    bin_size: float = 0.01,
                    outlier_rejection: OutlierRejectionMethod = OutlierRejectionMethod.SIGMA_CLIPPING,
                    weighting: WeightingScheme = WeightingScheme.EVEN,
                    normalization: AvgNormalization = AvgNormalization.NONE,
                    sigma_low: float = 3.0,
                    sigma_high: float = 3.0,
                    percentile: float = 0.9,
                    below_threshold_cutoff: float = 0.2,
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """Average multiple spectra with binning, outlier rejection, and weighting.

    Args:
        mz_arrays: list of m/z arrays
        intensity_arrays: list of intensity arrays
        bin_size: m/z bin width in Da (default 0.01)
        outlier_rejection: method for rejecting outlier intensities per bin
        weighting: per-spectrum weighting scheme
        normalization: intensity normalization before averaging
        sigma_low/sigma_high: sigma thresholds for clipping methods
        percentile: percentile for percentile clipping (default 0.9)
        below_threshold_cutoff: fraction cutoff for below_threshold rejection

    Returns:
        (averaged_mz, averaged_intensity) as numpy arrays
    """
    n_spectra = len(mz_arrays)
    if n_spectra == 0:
        return np.array([]), np.array([])
    if n_spectra == 1:
        return mz_arrays[0].copy(), intensity_arrays[0].copy()

    # Normalize
    norm_intensities, avg_tic = _normalize_spectra(intensity_arrays, normalization)

    # Compute weights
    weights = _compute_weights(intensity_arrays, weighting)

    # Find global m/z range
    all_mz = np.concatenate(mz_arrays)
    min_mz = all_mz.min()
    max_mz = all_mz.max()
    n_bins = int(np.ceil((max_mz - min_mz) / bin_size)) + 1

    # Bin peaks: for each bin, collect (intensity, spectrum_id, mz)
    bins_mz = [[] for _ in range(n_bins)]
    bins_intensity = [[] for _ in range(n_bins)]
    bins_spectra_id = [[] for _ in range(n_bins)]

    for spec_id in range(n_spectra):
        mz = mz_arrays[spec_id]
        inten = norm_intensities[spec_id]
        for j in range(len(mz)):
            idx = int((mz[j] - min_mz) / bin_size)
            if 0 <= idx < n_bins:
                bins_mz[idx].append(mz[j])
                bins_intensity[idx].append(inten[j])
                bins_spectra_id[idx].append(spec_id)

    # Zero-fill: for each bin, add zero entries for missing spectra
    for b in range(n_bins):
        if len(bins_spectra_id[b]) == 0:
            continue
        present = set(bins_spectra_id[b])
        if len(present) < n_spectra:
            avg_mz_in_bin = np.mean(bins_mz[b])
            for sid in range(n_spectra):
                if sid not in present:
                    bins_mz[b].append(avg_mz_in_bin)
                    bins_intensity[b].append(0.0)
                    bins_spectra_id[b].append(sid)

    # Get outlier rejection function
    reject_fn = _REJECTION_FUNCS.get(outlier_rejection, _reject_none)
    reject_kwargs = {
        'sigma_low': sigma_low, 'sigma_high': sigma_high,
        'percentile': percentile, 'cutoff': below_threshold_cutoff,
    }

    # Average each bin
    result_mz = []
    result_intensity = []

    for b in range(n_bins):
        if len(bins_intensity[b]) == 0:
            continue

        mz_vals = np.array(bins_mz[b])
        int_vals = np.array(bins_intensity[b])
        spec_ids = np.array(bins_spectra_id[b])

        # Outlier rejection on intensity values
        kept = reject_fn(int_vals, **reject_kwargs)
        if len(kept) == 0 or kept.sum() == 0:
            continue

        # Find which entries survived (match by value for simplicity)
        # For weighted average, we need spec_ids of surviving entries
        kept_set = set()
        temp_vals = int_vals.tolist()
        for v in kept:
            for i, tv in enumerate(temp_vals):
                if i not in kept_set and tv == v:
                    kept_set.add(i)
                    break

        kept_indices = sorted(kept_set)
        if not kept_indices:
            continue

        # Weighted average
        w_sum = 0.0
        int_sum = 0.0
        mz_sum = 0.0
        for i in kept_indices:
            w = weights[spec_ids[i]]
            int_sum += int_vals[i] * w
            mz_sum += mz_vals[i]
            w_sum += w

        if w_sum == 0 or int_sum == 0:
            continue

        result_mz.append(mz_sum / len(kept_indices))
        result_intensity.append(int_sum / w_sum)

    if not result_mz:
        return np.array([]), np.array([])

    out_mz = np.array(result_mz)
    out_int = np.array(result_intensity)

    # Restore absolute scale if normalized to TIC
    if normalization == AvgNormalization.ABSOLUTE_TIC:
        out_int *= avg_tic

    return out_mz, out_int
