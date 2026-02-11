# mzml-utils

Utilities for mzML mass spectrometry file processing, fragment ion calculation, and peak matching.

## Installation

```bash
pip install mzml-utils
```

For development:

```bash
git clone https://github.com/longpingfu/mzml-utils.git
cd mzml-utils
pip install -e ".[dev]"
```

## Modules

| Module | Description |
|--------|-------------|
| `reader` | mzML file I/O with indexed and sequential access |
| `ions` | Ion searching, mass matching, tolerance calculations |
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
