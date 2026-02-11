"""
Ion detection, mass matching, and tolerance calculations.

Provides functions for searching experimental spectra for target ions,
calculating mass errors, and ranking peak intensities.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Optional, Union


def ppm_error(observed: float, theoretical: float) -> float:
    """
    Calculate mass error in ppm.

    Args:
        observed: Observed m/z value.
        theoretical: Theoretical m/z value.

    Returns:
        Error in parts per million.
    """
    return (observed - theoretical) / theoretical * 1e6


def da_error(observed: float, theoretical: float) -> float:
    """
    Calculate mass error in Daltons.

    Args:
        observed: Observed m/z value.
        theoretical: Theoretical m/z value.

    Returns:
        Error in Da (observed - theoretical).
    """
    return observed - theoretical


def within_tolerance(observed: float, theoretical: float,
                     tolerance: float, unit: str = 'ppm') -> bool:
    """
    Check whether an observed m/z is within tolerance of a theoretical value.

    Args:
        observed: Observed m/z.
        theoretical: Theoretical m/z.
        tolerance: Tolerance value.
        unit: ``'ppm'`` or ``'Da'``.

    Returns:
        True if the absolute error is within the given tolerance.
    """
    if unit == 'ppm':
        return abs(ppm_error(observed, theoretical)) <= tolerance
    elif unit == 'Da':
        return abs(da_error(observed, theoretical)) <= tolerance
    else:
        raise ValueError(f"Unknown unit '{unit}'. Use 'ppm' or 'Da'.")


def find_ion(mz_array: np.ndarray,
             int_array: np.ndarray,
             target_mz: float,
             tolerance: float = 0.1,
             unit: str = 'Da') -> Optional[Dict]:
    """
    Find the most intense peak matching a target m/z within tolerance.

    Args:
        mz_array: Experimental m/z array.
        int_array: Experimental intensity array.
        target_mz: Target m/z to search for.
        tolerance: Mass tolerance.
        unit: ``'Da'`` or ``'ppm'``.

    Returns:
        Dict with ``mz``, ``intensity``, ``da_error``, ``ppm_error``
        for the best match, or None if no peak is found.
    """
    if len(mz_array) == 0:
        return None

    if unit == 'Da':
        mask = np.abs(mz_array - target_mz) <= tolerance
    elif unit == 'ppm':
        tol_da = target_mz * tolerance / 1e6
        mask = np.abs(mz_array - target_mz) <= tol_da
    else:
        raise ValueError(f"Unknown unit '{unit}'. Use 'ppm' or 'Da'.")

    if not mask.any():
        return None

    idx = np.argmax(int_array[mask])
    obs_mz = float(mz_array[mask][idx])
    obs_int = float(int_array[mask][idx])

    return {
        'mz': obs_mz,
        'intensity': obs_int,
        'da_error': da_error(obs_mz, target_mz),
        'ppm_error': ppm_error(obs_mz, target_mz),
    }


def search_ions(mz_array: np.ndarray,
                int_array: np.ndarray,
                ion_list: Dict[str, float],
                tolerance: float = 0.1,
                unit: str = 'Da',
                rel_threshold: float = 0.0) -> Dict[str, Optional[Dict]]:
    """
    Search a spectrum for multiple target ions.

    Args:
        mz_array: Experimental m/z array.
        int_array: Experimental intensity array.
        ion_list: Dict mapping ion names to target m/z values.
        tolerance: Mass tolerance.
        unit: ``'Da'`` or ``'ppm'``.
        rel_threshold: Minimum relative intensity (% of base peak)
            for a match to be reported. Default 0 (report all matches).

    Returns:
        Dict mapping each ion name to a match dict (with added ``rel_intensity``
        and ``above_threshold`` fields) or None if not found.
    """
    base_peak = float(int_array.max()) if len(int_array) > 0 else 0.0
    results: Dict[str, Optional[Dict]] = {}

    for name, target_mz in ion_list.items():
        match = find_ion(mz_array, int_array, target_mz, tolerance, unit)
        if match is not None:
            rel_int = match['intensity'] / base_peak * 100 if base_peak > 0 else 0.0
            match['rel_intensity'] = rel_int
            match['above_threshold'] = rel_int >= rel_threshold
            results[name] = match
        else:
            results[name] = None

    return results


def ion_rank(mz_array: np.ndarray,
             int_array: np.ndarray,
             target_mz: float,
             tolerance: float,
             unit: str = 'ppm') -> Optional[int]:
    """
    Find the intensity rank of an ion in the spectrum.

    Rank 1 = most intense peak overall.

    Args:
        mz_array: Experimental m/z array.
        int_array: Experimental intensity array.
        target_mz: Target m/z value.
        tolerance: Mass tolerance.
        unit: ``'ppm'`` or ``'Da'``.

    Returns:
        1-based intensity rank, or None if the ion is not found.
    """
    match = find_ion(mz_array, int_array, target_mz, tolerance, unit)
    if match is None:
        return None

    # Count how many peaks are more intense
    rank = int(np.sum(int_array > match['intensity'])) + 1
    return rank
