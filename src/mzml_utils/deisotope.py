"""
Averagine-scored deisotoping for centroided mass spectra.

Based on the MSFragger deisotoping algorithm (Teo et al., J. Proteome Res. 2021):
- Groups peaks into isotope clusters by m/z spacing
- Scores clusters against averagine distribution using KL divergence
- Validates subclusters to reject noise-inflated envelopes
- Retains singleton peaks (no isotope neighbors)

Uses Poisson approximation for averagine (no scipy dependency).
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .constants import (
    AVERAGINE_HEAVY_FRACTION,
    DEISOTOPE_SCORE_THRESHOLDS,
    NEUTRON_MASS,
    PROTON,
)


@dataclass
class IsotopeCluster:
    """An isotope cluster identified by deisotoping."""
    mono_idx: int              # index of monoisotopic peak in original arrays
    mono_mz: float             # monoisotopic m/z
    charge: int                # determined charge state
    n_peaks: int               # number of isotope peaks in the cluster
    peak_indices: List[int] = field(default_factory=list)  # indices of all peaks
    score: float = 0.0         # averagine KL score (0-1, higher = better)
    mono_intensity: float = 0.0  # intensity of monoisotopic peak


# ---------------------------------------------------------------------------
# Averagine isotope distribution (Poisson approximation)
# ---------------------------------------------------------------------------

def _poisson_distribution(neutral_mass: float, n_peaks: int) -> np.ndarray:
    """Theoretical isotope distribution using Poisson approximation of averagine.

    Parameters
    ----------
    neutral_mass : float
        Neutral mass of the ion in Da.
    n_peaks : int
        Number of isotope peaks to generate (M+0 through M+(n-1)).

    Returns
    -------
    np.ndarray
        Normalized probability distribution for isotope peaks.
    """
    lam = neutral_mass * AVERAGINE_HEAVY_FRACTION
    dist = np.array([
        math.exp(-lam + k * math.log(max(lam, 1e-10)) - math.lgamma(k + 1))
        for k in range(n_peaks)
    ])
    total = dist.sum()
    if total > 0:
        dist /= total
    return dist


def _kl_score(observed: np.ndarray, theoretical: np.ndarray) -> float:
    """Score isotope cluster against theoretical distribution.

    Uses KL divergence transformed to a 0-1 score: score = 1 / (1 + KL).
    Higher score = better match to averagine.

    Parameters
    ----------
    observed : np.ndarray
        Observed intensities (will be normalized).
    theoretical : np.ndarray
        Theoretical distribution (already normalized).

    Returns
    -------
    float
        Score in [0, 1]. 1.0 = perfect match.
    """
    n = min(len(observed), len(theoretical))
    obs = observed[:n].copy()
    theo = theoretical[:n].copy()

    # Normalize
    obs_sum = obs.sum()
    if obs_sum <= 0:
        return 0.0
    obs = obs / obs_sum

    # Avoid log(0) by adding small epsilon
    eps = 1e-10
    obs = np.clip(obs, eps, None)
    theo = np.clip(theo, eps, None)

    kl = np.sum(obs * np.log(obs / theo))
    return 1.0 / (1.0 + kl)


def _get_score_threshold(n_peaks: int,
                         thresholds: Dict[int, float]) -> float:
    """Get the minimum score threshold for a given cluster size."""
    if n_peaks in thresholds:
        return thresholds[n_peaks]
    # For clusters larger than max key, use the largest threshold
    max_key = max(thresholds.keys())
    if n_peaks > max_key:
        return thresholds[max_key]
    return 1.0  # shouldn't happen, but reject by default


# ---------------------------------------------------------------------------
# Core deisotoping
# ---------------------------------------------------------------------------

def deisotope(
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    max_charge: int = 4,
    tolerance_ppm: float = 15.0,
    min_score: Optional[Dict[int, float]] = None,
    return_clusters: bool = False,
) -> Union[
    Tuple[np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, List[IsotopeCluster]],
]:
    """Remove isotope peaks from a centroided spectrum, keeping monoisotopic peaks.

    Uses averagine-scored cluster detection with subcluster validation,
    inspired by the MSFragger deisotoping algorithm.

    Parameters
    ----------
    exp_mz : np.ndarray
        Experimental m/z values (must be sorted ascending).
    exp_intensity : np.ndarray
        Experimental intensity values.
    max_charge : int
        Maximum charge state to consider for isotope clusters.
    tolerance_ppm : float
        Tolerance in ppm for matching isotope peaks.
    min_score : dict, optional
        Minimum averagine KL scores by cluster size.
        Defaults to MSFragger thresholds: {2: 0.05, 3: 0.1, 4: 0.2, 5: 0.4, 6: 0.6}.
    return_clusters : bool
        If True, also return the list of identified IsotopeCluster objects.

    Returns
    -------
    filtered_mz : np.ndarray
        m/z values with isotope peaks removed.
    filtered_intensity : np.ndarray
        Intensity values with isotope peaks removed.
    clusters : list of IsotopeCluster, optional
        Only returned if return_clusters=True.
    """
    if min_score is None:
        min_score = DEISOTOPE_SCORE_THRESHOLDS

    if len(exp_mz) == 0:
        if return_clusters:
            return exp_mz, exp_intensity, []
        return exp_mz, exp_intensity

    n_peaks = len(exp_mz)
    is_isotope = np.zeros(n_peaks, dtype=bool)
    assigned = np.zeros(n_peaks, dtype=bool)  # assigned to a cluster (mono or isotope)
    clusters: List[IsotopeCluster] = []

    for i in range(n_peaks):
        if assigned[i]:
            continue

        best_cluster = None
        best_score = -1.0

        for z in range(1, max_charge + 1):
            spacing = NEUTRON_MASS / z

            # Collect consecutive isotope peaks
            cluster_indices = [i]
            cluster_intensities = [exp_intensity[i]]
            current_mz = exp_mz[i]

            for iso_offset in range(1, 10):  # up to M+9
                target_mz = exp_mz[i] + iso_offset * spacing
                tol = target_mz * tolerance_ppm / 1e6

                # Binary search for candidate peaks
                lo = np.searchsorted(exp_mz, target_mz - tol)
                hi = np.searchsorted(exp_mz, target_mz + tol, side='right')

                found = False
                best_j = -1
                best_err = float('inf')
                for j in range(lo, hi):
                    if j == i or is_isotope[j]:
                        continue
                    err = abs(exp_mz[j] - target_mz)
                    if err < best_err:
                        best_err = err
                        best_j = j

                if best_j >= 0:
                    cluster_indices.append(best_j)
                    cluster_intensities.append(exp_intensity[best_j])
                else:
                    break  # stop at first missing isotope

            if len(cluster_indices) < 2:
                continue

            # Sanity check: M+0 should not be negligible vs max peak
            # For peptides/glycopeptides (<5000 Da), monoisotopic is always
            # a significant fraction of the envelope
            max_int = max(cluster_intensities)
            if max_int > 0 and cluster_intensities[0] / max_int < 0.02:
                continue

            # Estimate neutral mass and score against averagine
            neutral_mass = (exp_mz[i] * z) - (z * PROTON)
            if neutral_mass <= 0:
                continue

            obs = np.array(cluster_intensities)
            theo = _poisson_distribution(neutral_mass, len(cluster_indices))
            score = _kl_score(obs, theo)

            threshold = _get_score_threshold(len(cluster_indices), min_score)
            if score < threshold:
                continue

            # Subcluster validation: all sub-envelopes [M+0..M+j] must pass
            passed_subcluster = True
            for sub_end in range(1, len(cluster_indices)):
                sub_obs = obs[:sub_end + 1]
                sub_theo = _poisson_distribution(neutral_mass, sub_end + 1)
                sub_score = _kl_score(sub_obs, sub_theo)
                sub_threshold = _get_score_threshold(sub_end + 1, min_score)
                if sub_end + 1 >= 2 and sub_score < sub_threshold:
                    # Trim cluster to the point before failure
                    cluster_indices = cluster_indices[:sub_end]
                    cluster_intensities = cluster_intensities[:sub_end]
                    obs = np.array(cluster_intensities)
                    if len(cluster_indices) < 2:
                        passed_subcluster = False
                    break

            if not passed_subcluster or len(cluster_indices) < 2:
                continue

            # Re-score trimmed cluster
            if len(cluster_indices) < len(cluster_intensities):
                obs = np.array(cluster_intensities[:len(cluster_indices)])
            theo = _poisson_distribution(neutral_mass, len(cluster_indices))
            score = _kl_score(obs, theo)

            if score > best_score:
                best_score = score
                best_cluster = IsotopeCluster(
                    mono_idx=i,
                    mono_mz=exp_mz[i],
                    charge=z,
                    n_peaks=len(cluster_indices),
                    peak_indices=list(cluster_indices),
                    score=score,
                    mono_intensity=exp_intensity[i],
                )

        if best_cluster is not None:
            clusters.append(best_cluster)
            assigned[best_cluster.peak_indices[0]] = True
            for idx in best_cluster.peak_indices[1:]:
                is_isotope[idx] = True
                assigned[idx] = True

    keep = ~is_isotope
    if return_clusters:
        return exp_mz[keep], exp_intensity[keep], clusters
    return exp_mz[keep], exp_intensity[keep]


# ---------------------------------------------------------------------------
# MS1 isotope envelope scoring
# ---------------------------------------------------------------------------

def score_ms1_envelope(
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    precursor_mz: float,
    charge: int,
    tolerance_ppm: float = 10.0,
    max_isotopes: int = 8,
) -> IsotopeCluster:
    """Score the MS1 isotope envelope around a precursor.

    Determines the true monoisotopic peak by walking down from the reported
    precursor m/z, then scores the full envelope against averagine.

    Parameters
    ----------
    exp_mz : np.ndarray
        MS1 m/z values (sorted ascending).
    exp_intensity : np.ndarray
        MS1 intensity values.
    precursor_mz : float
        Reported precursor m/z (may be M+0, M+1, M+2, etc.).
    charge : int
        Precursor charge state.
    tolerance_ppm : float
        Tolerance for peak matching.
    max_isotopes : int
        Maximum number of isotope peaks to collect.

    Returns
    -------
    IsotopeCluster
        Cluster with monoisotopic m/z, score, and peak indices.
        If no envelope found, returns cluster with score=0.
    """
    result = IsotopeCluster(
        mono_idx=-1,
        mono_mz=precursor_mz,
        charge=charge,
        n_peaks=0,
    )

    if len(exp_mz) == 0 or charge < 1:
        return result

    spacing = NEUTRON_MASS / charge

    def _find_peak(target_mz: float) -> Tuple[int, float, float]:
        """Find closest peak within tolerance. Returns (idx, mz, intensity)."""
        tol = target_mz * tolerance_ppm / 1e6
        lo = np.searchsorted(exp_mz, target_mz - tol)
        hi = np.searchsorted(exp_mz, target_mz + tol, side='right')
        if lo >= hi:
            return -1, 0.0, 0.0
        # Pick the closest peak
        segment = exp_mz[lo:hi]
        best_local = np.argmin(np.abs(segment - target_mz))
        idx = lo + best_local
        return idx, exp_mz[idx], exp_intensity[idx]

    # Step 1: Walk down from reported precursor to find monoisotopic
    mono_mz = precursor_mz
    mono_idx = -1
    mono_offset = 0  # how many isotope units below reported precursor

    # First find the reported precursor peak
    idx, mz, intensity = _find_peak(precursor_mz)
    if idx < 0:
        return result

    mono_mz = mz
    mono_idx = idx

    for offset in range(1, 6):
        candidate_mz = precursor_mz - offset * spacing
        c_idx, c_mz, c_int = _find_peak(candidate_mz)

        if c_idx < 0:
            # No peak found below — previous position is monoisotopic
            break

        # Check if peak below is significant (>5% of current best)
        ref_int = exp_intensity[mono_idx] if mono_idx >= 0 else intensity
        if c_int < ref_int * 0.05:
            break

        mono_mz = c_mz
        mono_idx = c_idx
        mono_offset = offset

    # Step 2: Walk up from monoisotopic to collect full envelope
    cluster_indices = [mono_idx]
    cluster_intensities = [exp_intensity[mono_idx]]

    for iso in range(1, max_isotopes):
        target_mz = mono_mz + iso * spacing
        idx, mz, intensity = _find_peak(target_mz)
        if idx < 0:
            break
        cluster_indices.append(idx)
        cluster_intensities.append(intensity)

    # Step 3: Score against averagine
    neutral_mass = (mono_mz * charge) - (charge * PROTON)
    if neutral_mass <= 0 or len(cluster_indices) < 2:
        result.mono_idx = mono_idx if mono_idx >= 0 else -1
        result.mono_mz = mono_mz
        result.n_peaks = len(cluster_indices)
        result.peak_indices = cluster_indices
        result.mono_intensity = cluster_intensities[0] if cluster_intensities else 0.0
        return result

    obs = np.array(cluster_intensities)
    theo = _poisson_distribution(neutral_mass, len(cluster_indices))
    score = _kl_score(obs, theo)

    result.mono_idx = mono_idx
    result.mono_mz = mono_mz
    result.charge = charge
    result.n_peaks = len(cluster_indices)
    result.peak_indices = cluster_indices
    result.score = score
    result.mono_intensity = cluster_intensities[0]

    return result
