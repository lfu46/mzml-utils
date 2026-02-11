"""
Scan pairing logic for MS1 cycle organization and HCD-EThcD pairing.

Groups MS2 scans into MS1 duty cycles and pairs triggered EThcD scans
with their corresponding HCD scans by precursor m/z matching.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .constants import NEUTRON_MASS


def match_precursors(mz1: float,
                     mz2: float,
                     charge: int,
                     tolerance_ppm: float = 10.0,
                     isotope_offsets: Optional[List[int]] = None
                     ) -> Tuple[bool, float, int]:
    """
    Check whether two precursor m/z values match, allowing isotope offsets.

    The instrument may select different isotope peaks for the same precursor
    in HCD vs EThcD scans. This function tries each isotope offset and
    returns the best match.

    Args:
        mz1: First precursor m/z (reference).
        mz2: Second precursor m/z (to be adjusted).
        charge: Precursor charge state.
        tolerance_ppm: Tolerance in ppm.
        isotope_offsets: Isotope offsets to try (default [0, 1, -1, 2, -2]).

    Returns:
        Tuple of (matched, ppm_error, best_isotope_offset).
    """
    if isotope_offsets is None:
        isotope_offsets = [0, 1, -1, 2, -2]

    best_ppm = float('inf')
    best_offset = 0

    for offset in isotope_offsets:
        adjusted = mz2 + offset * NEUTRON_MASS / charge if charge > 0 else mz2
        ppm = abs(mz1 - adjusted) / adjusted * 1e6
        if ppm < best_ppm:
            best_ppm = ppm
            best_offset = offset

    return best_ppm < tolerance_ppm, best_ppm, best_offset


def group_ms1_cycles(scan_list: List[dict]) -> List[Tuple[Optional[dict], List[dict]]]:
    """
    Group scans into MS1 duty cycles.

    Each cycle consists of one MS1 scan followed by its dependent MS2 scans.
    Scans must have an ``activation_type`` key (or ``act_type``).

    Args:
        scan_list: List of scan dicts, each with at least ``activation_type``
            (or ``act_type``) and ``scan_num``.

    Returns:
        List of (ms1_scan, [ms2_scans]) tuples. The first cycle may have
        ms1_scan=None if the file begins with MS2 scans.
    """
    cycles: List[Tuple[Optional[dict], List[dict]]] = []
    current_ms1: Optional[dict] = None
    current_ms2: List[dict] = []

    for s in scan_list:
        act = s.get('activation_type') or s.get('act_type', '')
        if act == 'MS1':
            if current_ms1 is not None or current_ms2:
                cycles.append((current_ms1, current_ms2))
            current_ms1 = s
            current_ms2 = []
        else:
            current_ms2.append(s)

    if current_ms1 is not None or current_ms2:
        cycles.append((current_ms1, current_ms2))

    return cycles


def pair_hcd_ethcd(cycles: List[Tuple[Optional[dict], List[dict]]],
                   tolerance_ppm: float = 10.0
                   ) -> Dict[int, int]:
    """
    Pair EThcD scans to their triggering HCD scans within MS1 cycles.

    Within each duty cycle, each EThcD scan is matched to the HCD scan
    with the closest precursor m/z (within tolerance, allowing isotope
    offsets).

    Args:
        cycles: Output of :func:`group_ms1_cycles`.
        tolerance_ppm: Precursor matching tolerance in ppm.

    Returns:
        Dict mapping ``hcd_scan_num -> ethcd_scan_num`` for paired scans.
    """
    paired: Dict[int, int] = {}

    for _ms1, ms2_scans in cycles:
        hcd_list = [s for s in ms2_scans
                    if (s.get('activation_type') or s.get('act_type')) == 'HCD']
        ethcd_list = [s for s in ms2_scans
                      if (s.get('activation_type') or s.get('act_type')) == 'EThcD']

        for ethcd in ethcd_list:
            best_match = None
            best_ppm = float('inf')

            ethcd_mz = ethcd.get('precursor_mz') or ethcd.get('prec_mz', 0)
            ethcd_charge = ethcd.get('precursor_charge') or ethcd.get('charge', 0)

            for hcd in hcd_list:
                hcd_mz = hcd.get('precursor_mz') or hcd.get('prec_mz', 0)
                matched, ppm, _offset = match_precursors(
                    ethcd_mz, hcd_mz, ethcd_charge, tolerance_ppm
                )
                if matched and ppm < best_ppm:
                    best_ppm = ppm
                    best_match = hcd

            if best_match is not None:
                hcd_scan = best_match.get('scan_num', 0)
                ethcd_scan = ethcd.get('scan_num', 0)
                paired[hcd_scan] = ethcd_scan

    return paired
