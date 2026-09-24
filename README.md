# mzml-utils

[![PyPI](https://img.shields.io/pypi/v/mzml-utils.svg)](https://pypi.org/project/mzml-utils/)
[![CI](https://github.com/lfu46/mzml-utils/actions/workflows/ci.yml/badge.svg)](https://github.com/lfu46/mzml-utils/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Utilities for mzML mass spectrometry file processing, fragment ion calculation, and peak matching.

## Installation

```bash
pip install mzml-utils
```

For development:

```bash
git clone https://github.com/lfu46/mzml-utils.git
cd mzml-utils
pip install -e ".[dev]"
```

## Modules

| Module | Description |
|--------|-------------|
| `reader` | mzML file I/O with indexed and sequential access |
| `ions` | Ion searching, mass matching, tolerance calculations |
| `xic` | Extracted ion chromatograms (intensity vs retention time) over a run; isotope co-elution evidence for a precursor; peak boundaries and peak-quality metrics |
| `chrom_extract` | Peak table → chromatograms across many runs, every row column carried through (the Bioconductor `chromExtract` analogue) |
| `isotopes` | Exact-composition isotope distributions, for envelopes the averagine cannot express (halogen tags, heavy labels) |
| `pairing` | MS1 duty-cycle grouping, HCD-EThcD scan pairing |
| `fragments` | Fragment ion calculator (b/y/c/z/Y/oxonium), peak matching, false match rate |
| `constants` | Physical constants, amino acid masses, common ion lists |
| `utils` | Spectrum ID parsing, file finding, modification string parsing |

## Quick Start

### Read an mzML file

```python
from mzml_utils import MzMLReader

with MzMLReader("experiment.mzML") as reader:
    spec = reader.get_spectrum(12345)
    print(f"Scan {spec.scan_num}: {spec.n_peaks} peaks, {spec.activation_type}")
    print(f"Precursor: {spec.precursor_mz:.4f} m/z, charge {spec.precursor_charge}")
```

### Search for diagnostic ions

```python
from mzml_utils import search_ions, OXONIUM_IONS

results = search_ions(spec.mz, spec.intensity, OXONIUM_IONS,
                      tolerance=0.1, unit='Da', rel_threshold=5.0)

for name, match in results.items():
    if match and match['above_threshold']:
        print(f"  {name}: {match['mz']:.4f} ({match['rel_intensity']:.1f}%)")
```

### Extract ion chromatograms (XIC)

```python
from mzml_utils import extract_xic, extract_xics, OXONIUM_IONS

# One precursor / (glyco)peptide XIC from MS1 (source = a path or an open reader)
xic = extract_xic("experiment.mzML", 1204.564, ms_level=1,
                  tolerance=20.0, unit='ppm', rt_range=(62.6, 74.1))
print(f"NL {xic.max_intensity:.2e} at RT {xic.apex_rt:.2f} min")

# Several oxonium ions + TIC/base-peak in a single pass (ChromatogramSet)
cs = extract_xics("experiment.mzML",
                  {k: OXONIUM_IONS[k] for k in ("HexNAc", "NeuAc")},
                  ms_level=2, activation='HCD', tolerance=20.0, unit='ppm')
print(cs.rt.shape, cs.tic.shape, cs["HexNAc"].max_intensity)
```

`ms_level=1` gives precursor XICs; `ms_level=2` gives oxonium/fragment XICs (add
`activation='HCD'` to skip ETD scans, `precursor_mz=...` for a single precursor).
RT is the reader's native unit (minutes for Thermo msconvert mzML). The `ms_level` and
`rt_range` filters are pushed into the reader (`iter_spectra(ms_level=, rt_range=)`), so a
spectrum cache never decodes the scans an MS1 XIC would drop.

### Check a precursor's isotopes at the trace level

```python
from pyteomics import mass
from mzml_utils import isotope_distribution, isotope_trace_evidence

# Do M+1 / M+2 co-elute with M, and does anything co-elute one isotope BELOW it?
ev = isotope_trace_evidence("experiment.mzML", 991.5001, 4, rt_center=75.9)
print(ev.coelution, ev.pattern_cosine, ev.m_minus_1_coelutes)

# An envelope the averagine cannot express: a 4-bromo-Phe peptide (M+2 > M)
comp = mass.Composition(sequence="TPENFPSK") + mass.Composition({"Br": 1, "H": -1})
dist = isotope_distribution(comp, n_peaks=3, charge=2)
print(dist.mz, dist.abundance)          # pass distribution=dist to isotope_trace_evidence
```

### Extract a peak table across many runs, with peak metrics

```python
import pandas as pd
from mzml_utils import extract_chromatograms, chromatogram_rows

# One row per target: m/z, RT window, run name, plus any annotation columns
table = psms[["psg", "run", "mz", "rt_min", "rt_max", "ms2_rt"]].to_dict("records")
res = extract_chromatograms(table, mzml_paths, by="run", apex_rt="ms2_rt",
                            tolerance=10.0, nearest_valley=True)
df = pd.DataFrame(chromatogram_rows(res))   # row columns + apex, boundaries, area, FWHM, ...
```

Each row gets its own m/z and RT window, and one pass per run serves every row of that
run. Without `by=`, every row is extracted from every run (cross-sample overlays).
`peak_metrics` / `peak_boundary` port the Chromatograms, MsQuality and MetaboCoreUtils R
code (boundaries, S/N, prominence, Kumler beta shape) and are pinned against those R
functions. `apex_rt=` measures the peak whose boundaries contain the MS2 time, and
`nearest_valley=True` stops R's left boundary from running to the start of a zero-filled
trace. These are measurements, not verdicts: no threshold is applied.

### Calculate fragment ions

```python
from mzml_utils import FragmentCalculator, match_peaks, parse_modifications

mods = parse_modifications("N-term(229.1629),4S(528.2859),19K(229.1629)")

calc = FragmentCalculator("AGYSQGATQYTQAQQTR", mods, precursor_charge=3)
theoretical = calc.get_all_ions_flat()

matched = match_peaks(theoretical, spec.mz, spec.intensity, tolerance_ppm=20.0)
print(f"Matched {len(matched)} of {len(theoretical)} theoretical ions")
```

### False match rate estimation

```python
from mzml_utils import calculate_false_match_rate

fmr = calculate_false_match_rate(theoretical, spec.mz, spec.intensity,
                                  tolerance_ppm=20.0)
print(f"FMR (peaks): {fmr.fmr_peaks*100:.1f}%")
print(f"FMR (intensity): {fmr.fmr_intensity*100:.1f}%")
```

### Pair HCD and EThcD scans

```python
from mzml_utils import MzMLReader, group_ms1_cycles, pair_hcd_ethcd

with MzMLReader("experiment.mzML") as reader:
    scans = [
        {
            'scan_num': s.scan_num,
            'activation_type': s.activation_type,
            'precursor_mz': s.precursor_mz,
            'precursor_charge': s.precursor_charge,
        }
        for s in reader.iter_spectra()
    ]

cycles = group_ms1_cycles(scans)
paired = pair_hcd_ethcd(cycles, tolerance_ppm=10.0)
# paired: {hcd_scan_num: ethcd_scan_num}
```

## Mass tolerance helpers

```python
from mzml_utils import ppm_error, da_error, within_tolerance

print(ppm_error(204.087, 204.0864))      # ~2.9 ppm
print(within_tolerance(204.087, 204.0864, 20, 'ppm'))  # True
```

## License

MIT
