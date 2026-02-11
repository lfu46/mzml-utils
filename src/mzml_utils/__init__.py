"""
mzml-utils: Utilities for mzML mass spectrometry file processing.

Modules:
    reader     - mzML file I/O (indexed + sequential)
    ions       - Ion search, mass matching, tolerance calculations
    pairing    - MS1 cycle organisation, HCD-EThcD pairing
    fragments  - Fragment ion calculator, peak matching, false match rate
    constants  - Physical constants, amino acid masses, common ion lists
    utils      - Spectrum ID parsing, file finding, modification parsing
"""

__version__ = "0.1.0"

# reader
from .reader import MzMLReader, Spectrum, classify_activation

# ions
from .ions import (
    ppm_error,
    da_error,
    within_tolerance,
    find_ion,
    search_ions,
    ion_rank,
)

# pairing
from .pairing import match_precursors, group_ms1_cycles, pair_hcd_ethcd

# fragments
from .fragments import (
    FragmentCalculator,
    TheoreticalIon,
    MatchedIon,
    FalseMatchRate,
    Modification,
    match_peaks,
    calculate_false_match_rate,
    calculate_annotation_statistics,
)

# constants (commonly used)
from .constants import (
    PROTON,
    H2O,
    NH3,
    NEUTRON_MASS,
    AA_MASSES,
    MOD_MASSES,
    OXONIUM_IONS,
)

# utils
from .utils import parse_spectrum_id, find_mzml_file, parse_modifications

__all__ = [
    # version
    "__version__",
    # reader
    "MzMLReader",
    "Spectrum",
    "classify_activation",
    # ions
    "ppm_error",
    "da_error",
    "within_tolerance",
    "find_ion",
    "search_ions",
    "ion_rank",
    # pairing
    "match_precursors",
    "group_ms1_cycles",
    "pair_hcd_ethcd",
    # fragments
    "FragmentCalculator",
    "TheoreticalIon",
    "MatchedIon",
    "FalseMatchRate",
    "Modification",
    "match_peaks",
    "calculate_false_match_rate",
    "calculate_annotation_statistics",
    # constants
    "PROTON",
    "H2O",
    "NH3",
    "NEUTRON_MASS",
    "AA_MASSES",
    "MOD_MASSES",
    "OXONIUM_IONS",
    # utils
    "parse_spectrum_id",
    "find_mzml_file",
    "parse_modifications",
]
