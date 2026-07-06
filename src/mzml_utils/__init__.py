"""
mzml-utils: Utilities for mzML mass spectrometry file processing.

Modules:
    reader     - mzML file I/O (indexed + sequential)
    ions       - Ion search, mass matching, tolerance calculations
    pairing    - MS1 cycle organisation, HCD-EThcD pairing
    fragments  - Fragment ion calculator, peak matching, false match rate
    deisotope  - Averagine-scored deisotoping and MS1 envelope scoring
    similarity - Spectral similarity metrics (cosine, entropy, KL, etc.)
    averaging  - Spectral averaging with outlier rejection and weighting
    proteases  - Protease definitions and in silico protein digestion
    constants  - Physical constants, amino acid masses, common ion lists
    utils      - Spectrum ID parsing, file finding, modification parsing
"""

__version__ = "0.1.0"

# reader
from .reader import MzMLReader, Spectrum, classify_activation

# spectrum cache -- cache-aware reader. open_spectra(path) is the SINGLE entry point
# for reading spectra: it returns a fast local SpectrumCache when a co-located
# <mzml_dir>/spectra_cache/<stem>.spectra.db exists, else a real MzMLReader (identical
# behaviour, same Spectrum objects). Build a cache once with build_cache().
from .spectrum_cache import (
    open_spectra,
    SpectrumCache,
    cache_path_for,
    build_cache,
    verify_cache,
)

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
from .pairing import (
    match_precursors,
    group_ms1_cycles,
    pair_hcd_ethcd,
    bin_and_sum_spectra,
    write_combined_etd_mzml,
)

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

# deisotope
from .deisotope import deisotope, score_ms1_envelope, IsotopeCluster

# canonical match layer
from .match_context import (
    MatchContext,
    canonical_match,
    from_spectrum,
    from_arrays,
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

# similarity
from .similarity import (
    spectral_similarity,
    pair_peaks,
    NormalizationScheme,
    SIMILARITY_METHODS,
)

# averaging
from .averaging import (
    average_spectra,
    OutlierRejectionMethod,
    WeightingScheme,
    AvgNormalization,
)

# proteases
from .proteases import (
    Protease,
    CleavageRule,
    dual_enzyme_digest,
    PROTEASES,
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
    # spectrum cache (cache-aware reader)
    "open_spectra",
    "SpectrumCache",
    "cache_path_for",
    "build_cache",
    "verify_cache",
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
    "bin_and_sum_spectra",
    "write_combined_etd_mzml",
    # fragments
    "FragmentCalculator",
    "TheoreticalIon",
    "MatchedIon",
    "FalseMatchRate",
    "Modification",
    "match_peaks",
    "calculate_false_match_rate",
    "calculate_annotation_statistics",
    # deisotope
    "deisotope",
    "score_ms1_envelope",
    "IsotopeCluster",
    # constants
    "PROTON",
    "H2O",
    "NH3",
    "NEUTRON_MASS",
    "AA_MASSES",
    "MOD_MASSES",
    "OXONIUM_IONS",
    # similarity
    "spectral_similarity",
    "pair_peaks",
    "NormalizationScheme",
    "SIMILARITY_METHODS",
    # averaging
    "average_spectra",
    "OutlierRejectionMethod",
    "WeightingScheme",
    "AvgNormalization",
    # proteases
    "Protease",
    "CleavageRule",
    "dual_enzyme_digest",
    "PROTEASES",
    # utils
    "parse_spectrum_id",
    "find_mzml_file",
    "parse_modifications",
]
