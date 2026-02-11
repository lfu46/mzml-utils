"""
Utility functions for spectrum ID parsing, file finding, and modification parsing.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple


def parse_spectrum_id(spectrum_str: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """
    Parse a spectrum identifier in FragPipe/MSFragger format.

    Format: ``{file_name}.{scan}.{scan}.{charge}``

    Args:
        spectrum_str: The spectrum identifier string.

    Returns:
        Tuple of (file_name, scan_number, charge).
        Returns (None, None, None) if parsing fails.
    """
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        try:
            return parts[0], int(parts[1]), int(parts[3])
        except (ValueError, IndexError):
            pass
    return None, None, None


def find_mzml_file(basename: str,
                   directory: str,
                   patterns: Optional[List[str]] = None) -> Optional[str]:
    """
    Find an mzML file in a directory, trying several naming conventions.

    Args:
        basename: Base file name (without extension or suffix).
        directory: Directory to search in.
        patterns: List of filename patterns to try. Each may contain
            ``{name}`` as a placeholder for *basename*. If None, uses
            common calibrated mzML patterns.

    Returns:
        Full path to the first matching file, or None.
    """
    if patterns is None:
        patterns = [
            "{name}_calibrated.mzML",
            "{name}_mz_calibrated.mzML",
            "{name}_ppm_calibrated.mzML",
            "{name}.mzML",
        ]

    for pattern in patterns:
        path = os.path.join(directory, pattern.format(name=basename))
        if os.path.exists(path):
            return path

    # Fallback: search for any file starting with basename
    try:
        for f in os.listdir(directory):
            if f.startswith(basename) and f.endswith('.mzML'):
                return os.path.join(directory, f)
    except OSError:
        pass

    return None


def parse_modifications(mod_string: str) -> List[Dict]:
    """
    Parse a modification string from FragPipe/OGlycoTM format.

    Format: ``"N-term(229.1629),4S(528.2859),19K(229.1629)"``

    Each modification is returned as a dict with:
      - ``position``: 0 for N-term, -1 for C-term, 1-based for residues
      - ``residue``: ``'N-term'``, ``'C-term'``, or amino acid letter
      - ``mass``: Modification mass as float

    Args:
        mod_string: Comma-separated modification string.

    Returns:
        List of modification dicts.
    """
    if not mod_string or mod_string == 'nan' or str(mod_string) == 'nan':
        return []

    mods = []
    for part in str(mod_string).split(','):
        part = part.strip()
        if not part or '(' not in part:
            continue

        try:
            mass_str = part[part.find('(') + 1:part.find(')')]
            mass = float(mass_str)
            position_part = part[:part.find('(')]

            if position_part == 'N-term':
                mods.append({'position': 0, 'residue': 'N-term', 'mass': mass})
            elif position_part == 'C-term':
                mods.append({'position': -1, 'residue': 'C-term', 'mass': mass})
            else:
                pos = ''
                res = ''
                for char in position_part:
                    if char.isdigit():
                        pos += char
                    else:
                        res += char
                if pos and res:
                    mods.append({'position': int(pos), 'residue': res, 'mass': mass})
        except (ValueError, IndexError):
            continue

    return mods
