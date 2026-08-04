"""
Spectral similarity metrics for comparing two mass spectra.

Implements 9 similarity measures based on mzLib (Smith lab, UW-Madison):
cosine, spectral contrast angle, euclidean, bray-curtis, pearson,
dot product, spectral entropy, KL divergence, and Searle similarity.

Reference: mzLib SpectralSimilarity.cs
"""

import numpy as np
import math
from enum import Enum
from typing import Optional, Tuple, List


class NormalizationScheme(Enum):
    MOST_ABUNDANT = "most_abundant"
    SPECTRUM_SUM = "spectrum_sum"
    SQRT_SPECTRUM_SUM = "sqrt_spectrum_sum"
    UNNORMALIZED = "unnormalized"


def _normalize(intensity: np.ndarray, scheme: NormalizationScheme) -> np.ndarray:
    """Normalize intensity array."""
    if scheme == NormalizationScheme.MOST_ABUNDANT:
        mx = intensity.max()
        return intensity / mx if mx > 0 else intensity.copy()
    elif scheme == NormalizationScheme.SPECTRUM_SUM:
        s = intensity.sum()
        return intensity / s if s > 0 else intensity.copy()
    elif scheme == NormalizationScheme.SQRT_SPECTRUM_SUM:
        sq = np.sqrt(intensity)
        s = sq.sum()
        return sq / s if s > 0 else sq
    else:
        return intensity.copy()


def pair_peaks(mz1: np.ndarray, intensity1: np.ndarray,
               mz2: np.ndarray, intensity2: np.ndarray,
               tolerance_ppm: float = 20.0,
               all_peaks: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Pair peaks between two spectra using greedy intensity-descending matching.

    Each peak can only be paired once. Higher-intensity peaks are matched first.

    Args:
        mz1, intensity1: first spectrum (e.g., experimental)
        mz2, intensity2: second spectrum (e.g., theoretical/library)
        tolerance_ppm: matching tolerance in ppm
        all_peaks: if True, include unpaired peaks from spectrum 1 as (intensity, 0)

    Returns:
        (paired_1, paired_2): aligned intensity arrays
    """
    n1, n2 = len(mz1), len(mz2)
    if n1 == 0 or n2 == 0:
        return np.array([]), np.array([])

    # Sort indices by intensity descending
    order1 = np.argsort(-intensity1)
    order2 = np.argsort(-intensity2)

    used1 = np.zeros(n1, dtype=bool)
    pairs_1 = []
    pairs_2 = []

    # For each peak in spectrum 2 (descending intensity), find best match in spectrum 1
    for j in order2:
        best_i = -1
        best_inten = -1
        for i in order1:
            if used1[i]:
                continue
            ppm = abs(mz1[i] - mz2[j]) / max(mz1[i], mz2[j]) * 1e6
            if ppm <= tolerance_ppm:
                if intensity1[i] > best_inten:
                    best_i = i
                    best_inten = intensity1[i]
        if best_i >= 0:
            pairs_1.append(intensity1[best_i])
            pairs_2.append(intensity2[j])
            used1[best_i] = True
        else:
            pairs_1.append(0.0)
            pairs_2.append(intensity2[j])

    # Add unpaired peaks from spectrum 1
    if all_peaks:
        for i in range(n1):
            if not used1[i]:
                pairs_1.append(intensity1[i])
                pairs_2.append(0.0)

    return np.array(pairs_1), np.array(pairs_2)


def spectral_similarity(mz1: np.ndarray, intensity1: np.ndarray,
                         mz2: np.ndarray, intensity2: np.ndarray,
                         method: str = "cosine",
                         tolerance_ppm: float = 20.0,
                         normalization: NormalizationScheme = NormalizationScheme.MOST_ABUNDANT,
                         min_mz: float = 0.0,
                         all_peaks: bool = True) -> Optional[float]:
    """Compute spectral similarity between two spectra.

    Args:
        mz1, intensity1: first spectrum
        mz2, intensity2: second spectrum
        method: one of 'cosine', 'spectral_contrast_angle', 'euclidean',
                'bray_curtis', 'pearson', 'dot_product', 'spectral_entropy',
                'kl_divergence', 'searle'
        tolerance_ppm: peak matching tolerance
        normalization: intensity normalization scheme
        min_mz: filter out peaks below this m/z
        all_peaks: include unpaired peaks in comparison

    Returns:
        Similarity score (higher = more similar, except KL divergence where lower = more similar).
        Returns None if no valid pairs.
    """
    # Filter by min_mz
    if min_mz > 0:
        mask1 = mz1 >= min_mz
        mz1, intensity1 = mz1[mask1], intensity1[mask1]
        mask2 = mz2 >= min_mz
        mz2, intensity2 = mz2[mask2], intensity2[mask2]

    # Filter zero intensity
    mask1 = intensity1 > 0
    mz1, intensity1 = mz1[mask1], intensity1[mask1]
    mask2 = intensity2 > 0
    mz2, intensity2 = mz2[mask2], intensity2[mask2]

    if len(mz1) == 0 or len(mz2) == 0:
        return None

    # Normalize
    norm1 = _normalize(intensity1, normalization)
    norm2 = _normalize(intensity2, normalization)

    # Pair peaks
    p, q = pair_peaks(mz1, norm1, mz2, norm2, tolerance_ppm, all_peaks)

    if len(p) == 0:
        return None

    # Compute metric
    return _compute_metric(p, q, method, norm1, norm2)


def _compute_metric(p: np.ndarray, q: np.ndarray, method: str,
                    full_1: np.ndarray = None, full_2: np.ndarray = None) -> Optional[float]:
    """Compute similarity metric from paired intensity arrays."""
    if len(p) == 0:
        return None

    method = method.lower().replace('-', '_').replace(' ', '_')

    if method == "cosine":
        num = np.dot(p, q)
        denom = np.sqrt(np.dot(p, p) * np.dot(q, q))
        return float(num / denom) if denom > 0 else 0.0

    elif method == "spectral_contrast_angle":
        cos = _compute_metric(p, q, "cosine")
        if cos is None:
            return None
        # Clamp to [-1, 1] for numerical stability
        cos = max(-1.0, min(1.0, cos))
        return float(1.0 - 2.0 * math.acos(cos) / math.pi)

    elif method == "euclidean":
        return float(1.0 - np.sqrt(np.sum((p - q) ** 2)))

    elif method == "bray_curtis":
        denom = np.sum(p + q)
        if denom == 0:
            return None
        return float(1.0 - np.sum(np.abs(p - q)) / denom)

    elif method == "pearson":
        mp, mq = p.mean(), q.mean()
        dp, dq = p - mp, q - mq
        denom = np.sqrt(np.dot(dp, dp) * np.dot(dq, dq))
        if denom == 0:
            return -1.0
        return float(np.dot(dp, dq) / denom)

    elif method == "dot_product":
        return float(np.dot(p, q))

    elif method == "spectral_entropy":
        # Requires spectrum_sum normalization for proper entropy calculation
        def entropy(x):
            x = x[x > 0]
            return float(-np.sum(x * np.log(x)))

        h_p = entropy(p[p > 0]) if np.any(p > 0) else 0.0
        h_q = entropy(q[q > 0]) if np.any(q > 0) else 0.0
        combined = (p + q) / 2.0
        h_combined = entropy(combined[combined > 0]) if np.any(combined > 0) else 0.0
        return float(1.0 - (2.0 * h_combined - h_p - h_q) / math.log(4))

    elif method == "kl_divergence":
        # Lower = more similar. Adds epsilon to avoid log(0).
        eps = 1e-9
        p_smooth = (p + eps) / (p + eps).sum()
        q_smooth = (q + eps) / (q + eps).sum()
        return float(np.sum(p_smooth * np.log(p_smooth / q_smooth)))

    elif method == "searle":
        # From Brian Searle, ASMS 2022. Best with sqrt_spectrum_sum normalization.
        diff_sq = np.sum((p - q) ** 2)
        if diff_sq == 0:
            return float('inf')  # perfect match
        return float(math.log(1.0 / diff_sq))

    else:
        raise ValueError(f"Unknown similarity method: {method}. Choose from: "
                         "cosine, spectral_contrast_angle, euclidean, bray_curtis, "
                         "pearson, dot_product, spectral_entropy, kl_divergence, searle")


# Convenience aliases
SIMILARITY_METHODS = [
    "cosine", "spectral_contrast_angle", "euclidean", "bray_curtis",
    "pearson", "dot_product", "spectral_entropy", "kl_divergence", "searle"
]
