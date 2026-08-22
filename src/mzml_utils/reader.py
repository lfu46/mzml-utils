"""
mzML file reading with indexed and sequential access.

Provides a clean abstraction over pyteomics mzML parsing with
scan metadata extraction and activation type classification.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Generator, List, Optional

from pyteomics import mzml


@dataclass
class Spectrum:
    """Container for a single mass spectrum."""
    scan_num: int
    mz: np.ndarray
    intensity: np.ndarray
    ms_level: int = 0
    rt: float = 0.0
    filter_string: str = ""
    precursor_mz: float = 0.0
    precursor_charge: int = 0
    precursor_intensity: float = 0.0
    activation_type: str = ""
    tic: float = 0.0
    base_peak_mz: float = 0.0
    base_peak_intensity: float = 0.0
    isolation_window_target: float = 0.0
    isolation_window_lower: float = 0.0
    isolation_window_upper: float = 0.0
    ion_injection_time: float = 0.0
    """Trap fill time in ms; 0.0 when the instrument did not report it.

    Instruments cap this **per scan type** -- on an Orbitrap Eclipse, MS1, HCD and
    EThcD each carry their own maximum -- so a figure pooled across scan types is
    meaningless. Group by ``activation_type`` before summarising.
    """

    @property
    def n_peaks(self) -> int:
        return len(self.mz)

    @property
    def isolation_width(self) -> float:
        """Total isolation window width in Da."""
        return self.isolation_window_lower + self.isolation_window_upper


def classify_activation(filter_string: str) -> str:
    """
    Classify the activation type from an instrument filter string.

    Args:
        filter_string: Raw filter string from the mzML scan metadata.

    Returns:
        One of 'MS1', 'HCD', 'EThcD', 'ETD', 'CID', or 'other'.
    """
    fs = filter_string.lower()
    if 'ms ' in fs and 'ms2' not in fs and '@' not in fs:
        return 'MS1'
    if '@etd' in fs and '@hcd' in fs:
        return 'EThcD'
    if '@etd' in fs:
        return 'ETD'
    if '@hcd' in fs:
        return 'HCD'
    if '@cid' in fs:
        return 'CID'
    return 'other'


def _parse_scan_num(spec: dict) -> int:
    """Extract scan number from a pyteomics spectrum dict."""
    scan_id = spec.get('id', '')
    if 'scan=' in scan_id:
        return int(scan_id.split('scan=')[-1])
    return 0


def _extract_spectrum(spec: dict) -> Spectrum:
    """Convert a pyteomics spectrum dict to a Spectrum object."""
    scan_num = _parse_scan_num(spec)
    ms_level = spec.get('ms level', 0)

    mz_array = spec.get('m/z array', np.array([]))
    int_array = spec.get('intensity array', np.array([]))

    rt = 0.0
    filter_string = ''
    inj_time = 0.0
    if 'scanList' in spec and 'scan' in spec['scanList']:
        scan_info = spec['scanList']['scan'][0]
        rt_val = scan_info.get('scan start time')
        if rt_val is not None:
            rt = float(rt_val)
        filter_string = scan_info.get('filter string', '')
        inj_val = scan_info.get('ion injection time')
        if inj_val is not None:
            inj_time = float(inj_val)

    prec_mz = 0.0
    prec_charge = 0
    prec_intensity = 0.0
    iso_target = 0.0
    iso_lower = 0.0
    iso_upper = 0.0
    if 'precursorList' in spec:
        prec = spec['precursorList']['precursor'][0]
        if 'selectedIonList' in prec:
            sel_ion = prec['selectedIonList']['selectedIon'][0]
            prec_mz = float(sel_ion.get('selected ion m/z', 0))
            prec_charge = int(sel_ion.get('charge state', 0))
            prec_intensity = float(sel_ion.get('peak intensity', 0))
        if 'isolationWindow' in prec:
            iso_win = prec['isolationWindow']
            iso_target = float(iso_win.get('isolation window target m/z', 0))
            iso_lower = float(iso_win.get('isolation window lower offset', 0))
            iso_upper = float(iso_win.get('isolation window upper offset', 0))

    if ms_level == 1:
        act_type = 'MS1'
    else:
        act_type = classify_activation(filter_string)

    return Spectrum(
        scan_num=scan_num,
        mz=mz_array,
        intensity=int_array,
        ms_level=ms_level,
        rt=rt,
        filter_string=filter_string,
        precursor_mz=prec_mz,
        precursor_charge=prec_charge,
        precursor_intensity=prec_intensity,
        activation_type=act_type,
        tic=float(spec.get('total ion current', 0)),
        base_peak_mz=float(spec.get('base peak m/z', 0)),
        base_peak_intensity=float(spec.get('base peak intensity', 0)),
        isolation_window_target=iso_target,
        isolation_window_lower=iso_lower,
        isolation_window_upper=iso_upper,
        ion_injection_time=inj_time,
    )


class MzMLReader:
    """
    Read mzML files with indexed or sequential access.

    Args:
        path: Path to the mzML file.
        use_index: Whether to use indexed access (default True).
            Set to False for sequential-only reading of non-indexed files.

    Example::

        reader = MzMLReader("experiment.mzML")
        spec = reader.get_spectrum(1234)
        print(spec.mz, spec.intensity, spec.precursor_mz)
        reader.close()

        # Or as a context manager:
        with MzMLReader("experiment.mzML") as reader:
            for spec in reader.iter_spectra():
                print(spec.scan_num, spec.activation_type)
    """

    def __init__(self, path: str, use_index: bool = True):
        self.path = path
        self._use_index = use_index
        self._reader = mzml.MzML(path, use_index=use_index)

    def close(self):
        """Close the underlying file reader."""
        if hasattr(self._reader, 'close'):
            self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def get_spectrum(self, scan_num: int) -> Optional[Spectrum]:
        """
        Retrieve a single spectrum by scan number (indexed access).

        Args:
            scan_num: The scan number to retrieve.

        Returns:
            A Spectrum object, or None if the scan is not found.
        """
        scan_id = f'controllerType=0 controllerNumber=1 scan={scan_num}'
        try:
            spec = self._reader.get_by_id(scan_id)
        except (KeyError, IndexError):
            return None
        return _extract_spectrum(spec)

    def iter_spectra(self) -> Generator[Spectrum, None, None]:
        """Iterate over all spectra sequentially."""
        with mzml.MzML(self.path) as reader:
            for spec in reader:
                yield _extract_spectrum(spec)

    def find_best_ms1(self, ms2_scan: int, precursor_mz: float,
                      n_ms1: int = 3, tol_da: float = 0.02) -> Optional[int]:
        """
        Find the MS1 scan with the strongest precursor peak.

        Searches backwards from *ms2_scan* through the preceding *n_ms1*
        MS1 scans and returns the one whose peak closest to *precursor_mz*
        has the highest intensity.

        Args:
            ms2_scan: Scan number of the MS2 spectrum.
            precursor_mz: Expected precursor m/z.
            n_ms1: Number of preceding MS1 scans to check (default 3).
            tol_da: Tolerance in Da for matching the precursor peak.

        Returns:
            Scan number of the best MS1, or the closest MS1 if the
            precursor is not detected in any of them.
        """
        ms1_scans: List[int] = []
        scan = ms2_scan - 1
        while len(ms1_scans) < n_ms1 and scan > 0:
            spec = self.get_spectrum(scan)
            if spec is not None and spec.ms_level == 1:
                ms1_scans.append(scan)
            scan -= 1

        best_scan = None
        best_int = 0.0
        for ms1_scan in ms1_scans:
            spec = self.get_spectrum(ms1_scan)
            if spec is None:
                continue
            mask = np.abs(spec.mz - precursor_mz) < tol_da
            if mask.any():
                max_int = float(spec.intensity[mask].max())
                if max_int > best_int:
                    best_int = max_int
                    best_scan = ms1_scan

        # Fallback to closest MS1 if precursor not found
        if best_scan is None and ms1_scans:
            best_scan = ms1_scans[0]
        return best_scan

    def extract_metadata(self, scan_num: int) -> Optional[Dict]:
        """
        Extract metadata for a single scan as a plain dict.

        Returns:
            Dict with keys: scan_num, ms_level, rt, filter_string,
            precursor_mz, precursor_charge, activation_type,
            isolation_window_target, isolation_window_lower,
            isolation_window_upper, isolation_width, ion_injection_time.
            None if the scan is not found.
        """
        spec = self.get_spectrum(scan_num)
        if spec is None:
            return None
        return {
            'scan_num': spec.scan_num,
            'ms_level': spec.ms_level,
            'rt': spec.rt,
            'filter_string': spec.filter_string,
            'precursor_mz': spec.precursor_mz,
            'precursor_charge': spec.precursor_charge,
            'activation_type': spec.activation_type,
            'isolation_window_target': spec.isolation_window_target,
            'isolation_window_lower': spec.isolation_window_lower,
            'isolation_window_upper': spec.isolation_window_upper,
            'isolation_width': spec.isolation_width,
            'ion_injection_time': spec.ion_injection_time,
        }
