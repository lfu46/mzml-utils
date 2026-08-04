"""
Scan pairing logic for MS1 cycle organization and HCD-EThcD pairing.

Groups MS2 scans into MS1 duty cycles and pairs triggered EThcD scans
with their corresponding HCD scans by precursor m/z matching.

Also provides combined ETD mzML generation: sums multiple triggered
ETD scans (e.g., SA 10/15/20%) into a single composite scan for
improved OPair localization.
"""

from __future__ import annotations

import numpy as np
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


def bin_and_sum_spectra(mz_arrays: List[np.ndarray],
                        intensity_arrays: List[np.ndarray],
                        ppm_tol: float = 20.0
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Combine multiple spectra by m/z binning with intensity summation.

    For peaks at the same m/z (within *ppm_tol*), keeps the m/z from
    the most intense observation and sums the intensities.

    Args:
        mz_arrays: List of m/z arrays to combine.
        intensity_arrays: List of intensity arrays to combine.
        ppm_tol: Tolerance in ppm for grouping peaks.

    Returns:
        Tuple of (combined_mz, combined_intensity) arrays.
    """
    if not mz_arrays:
        return np.array([]), np.array([])

    all_mz = np.concatenate(mz_arrays)
    all_int = np.concatenate(intensity_arrays)

    if len(all_mz) == 0:
        return np.array([]), np.array([])

    order = np.argsort(all_mz)
    all_mz = all_mz[order]
    all_int = all_int[order]

    binned_mz: List[float] = []
    binned_int: List[float] = []
    i = 0
    while i < len(all_mz):
        cluster_mz = [all_mz[i]]
        cluster_int = [all_int[i]]
        j = i + 1
        while j < len(all_mz):
            ppm = (all_mz[j] - cluster_mz[0]) / cluster_mz[0] * 1e6
            if ppm <= ppm_tol:
                cluster_mz.append(all_mz[j])
                cluster_int.append(all_int[j])
                j += 1
            else:
                break
        best_idx = int(np.argmax(cluster_int))
        binned_mz.append(cluster_mz[best_idx])
        binned_int.append(sum(cluster_int))
        i = j

    return np.array(binned_mz), np.array(binned_int)


def write_combined_etd_mzml(input_path: str,
                             output_path: str,
                             ppm_tol: float = 20.0) -> Dict:
    """Create an mzML with ETD triplets combined into single scans.

    For each HCD-triggered ETD triplet (e.g., SA 10/15/20%), sums the
    three ETD spectra into one combined scan. The first ETD scan keeps
    its original scan number with the combined spectrum; the 2nd and
    3rd scans are removed.

    Args:
        input_path: Path to the original mzML file.
        output_path: Path to write the combined mzML.
        ppm_tol: Tolerance in ppm for m/z binning during combination.

    Returns:
        Dict with keys: n_triplets, n_input_scans, n_output_scans.

    Requires ``psims`` for mzML writing::

        pip install psims

    Example::

        from mzml_utils import write_combined_etd_mzml

        stats = write_combined_etd_mzml(
            "experiment.mzML",
            "experiment_combined_ETD.mzML"
        )
        print(f"Combined {stats['n_triplets']} ETD triplets")
    """
    from pyteomics import mzml
    from psims.mzml import MzMLWriter

    # Read all spectra
    all_spectra = []
    with mzml.MzML(input_path) as reader:
        for spec in reader:
            all_spectra.append(spec)

    # Find ETD triplets
    triplets = []
    i = 0
    while i < len(all_spectra):
        spec = all_spectra[i]
        ms_level = spec.get('ms level', 0)
        fs = ''
        if 'scanList' in spec and 'scan' in spec['scanList']:
            fs = spec['scanList']['scan'][0].get('filter string', '')

        if ms_level == 2 and '@hcd' in fs.lower() and '@etd' not in fs.lower():
            prec_mz = 0.0
            if 'precursorList' in spec:
                sel = spec['precursorList']['precursor'][0].get('selectedIonList', {})
                ions = sel.get('selectedIon', [{}])
                prec_mz = float(ions[0].get('selected ion m/z', 0))

            etd_indices = []
            j = i + 1
            while j < len(all_spectra) and len(etd_indices) < 3:
                ns = all_spectra[j]
                nfs = ''
                if 'scanList' in ns and 'scan' in ns['scanList']:
                    nfs = ns['scanList']['scan'][0].get('filter string', '')
                nl = ns.get('ms level', 0)
                if nl == 1:
                    break
                if nl == 2 and 'etd' in nfs.lower():
                    n_mz = 0.0
                    if 'precursorList' in ns:
                        sel2 = ns['precursorList']['precursor'][0].get('selectedIonList', {})
                        ions2 = sel2.get('selectedIon', [{}])
                        n_mz = float(ions2[0].get('selected ion m/z', 0))
                    if prec_mz > 0 and abs(n_mz - prec_mz) < 0.5:
                        etd_indices.append(j)
                j += 1

            if len(etd_indices) == 3:
                triplets.append((i, etd_indices))
        i += 1

    # Build skip/combine maps
    skip_indices = set()
    combine_map = {}
    for _, etd_indices in triplets:
        combine_map[etd_indices[0]] = etd_indices
        skip_indices.add(etd_indices[1])
        skip_indices.add(etd_indices[2])

    n_output = len(all_spectra) - len(skip_indices)

    # Write combined mzML
    with MzMLWriter(open(output_path, 'wb'), close=True) as writer:
        writer.controlled_vocabularies()
        writer.file_description(
            file_contents=['MSn spectrum', 'centroid spectrum'],
            source_files=[])
        writer.software_list([
            {'id': 'combine_etd', 'version': '1.0',
             'params': ['custom unreleased software tool']}])
        writer.instrument_configuration_list([
            writer.InstrumentConfiguration(
                id='IC1',
                component_list=[
                    writer.Source(order=1, params=['electrospray ionization']),
                    writer.Analyzer(order=2, params=['orbitrap']),
                    writer.Detector(order=3, params=['inductive detector']),
                ])])
        writer.data_processing_list([{
            'id': 'combined_etd_processing',
            'processing_methods': [{
                'order': 1, 'software_reference': 'combine_etd',
                'params': ['peak picking']}]}])

        with writer.run(id='run1', instrument_configuration='IC1'):
            with writer.spectrum_list(count=n_output):
                for i, spec in enumerate(all_spectra):
                    if i in skip_indices:
                        continue

                    ms_level = spec.get('ms level', 1)
                    mz_array = spec.get('m/z array', np.array([]))
                    int_array = spec.get('intensity array', np.array([]))

                    if i in combine_map:
                        mz_arrays = [all_spectra[idx].get('m/z array', np.array([]))
                                     for idx in combine_map[i]]
                        int_arrays = [all_spectra[idx].get('intensity array', np.array([]))
                                      for idx in combine_map[i]]
                        mz_array, int_array = bin_and_sum_spectra(
                            mz_arrays, int_arrays, ppm_tol=ppm_tol)

                    scan_id = spec.get('id', '')
                    scan_num = int(scan_id.split('scan=')[-1]) if 'scan=' in scan_id else i

                    rt = 0.0
                    filter_string = ''
                    if 'scanList' in spec and 'scan' in spec['scanList']:
                        rt_val = spec['scanList']['scan'][0].get('scan start time')
                        if rt_val is not None:
                            rt = float(rt_val)
                        filter_string = spec['scanList']['scan'][0].get('filter string', '')

                    spec_params = [{'ms level': ms_level}, 'centroid spectrum', 'positive scan']
                    if spec.get('total ion current'):
                        spec_params.append({'total ion current': float(spec['total ion current'])})
                    if spec.get('base peak m/z'):
                        spec_params.append({'base peak m/z': float(spec['base peak m/z'])})
                    if spec.get('base peak intensity'):
                        spec_params.append({'base peak intensity': float(spec['base peak intensity'])})

                    precursor_info = None
                    if ms_level >= 2 and 'precursorList' in spec:
                        prec = spec['precursorList']['precursor'][0]
                        sel_ions = prec.get('selectedIonList', {}).get('selectedIon', [{}])
                        if sel_ions:
                            sel = sel_ions[0]
                            p_mz = float(sel.get('selected ion m/z', 0))
                            p_z = int(sel.get('charge state', 0))
                            p_int = float(sel.get('peak intensity', 0))
                            iso_win = prec.get('isolationWindow', {})
                            iso_t = float(iso_win.get('isolation window target m/z', p_mz))
                            iso_lo = float(iso_win.get('isolation window lower offset', 0.7))
                            iso_hi = float(iso_win.get('isolation window upper offset', 0.7))
                            act = prec.get('activation', {})
                            act_params = []
                            if 'beam-type collision-induced dissociation' in act:
                                act_params.append('beam-type collision-induced dissociation')
                            if 'electron transfer dissociation' in act:
                                act_params.append('electron transfer dissociation')
                            if 'collision energy' in act:
                                act_params.append({'collision energy': float(act['collision energy'])})
                            precursor_info = {
                                'mz': p_mz, 'intensity': p_int, 'charge': p_z,
                                'isolation_window_args': {
                                    'target': iso_t, 'lower': iso_lo, 'upper': iso_hi},
                                'activation': act_params or ['beam-type collision-induced dissociation'],
                            }

                    writer.write_spectrum(
                        mz_array=mz_array, intensity_array=int_array,
                        id=f"controllerType=0 controllerNumber=1 scan={scan_num}",
                        params=spec_params, scan_start_time=rt,
                        precursor_information=precursor_info,
                        scan_params=[{'filter string': filter_string}] if filter_string else [])

    return {
        'n_triplets': len(triplets),
        'n_input_scans': len(all_spectra),
        'n_output_scans': n_output,
    }
