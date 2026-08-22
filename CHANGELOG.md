# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `Spectrum.ion_injection_time`, populated by both the indexed mzML reader and the spectrum cache (#2).
- `XIC.fwhm` and `XIC.points_across_peak` for chromatographic peak-width QC (#3).
- Isobaric reporter-ion extraction (`mzml_utils.reporters`).
- `mzml_utils.structure`: one resolver for every structure database (AlphaFold DB, PDBe,
  SWISS-MODEL, AlphaFill) via 3D-Beacons, plus UniProt features, InterPro domains and
  UniProt↔PDB residue mapping.
- `mzml_utils.open_spectra()` — the single cache-aware reader: a local `SpectrumCache` when a
  co-located `spectra_cache/<stem>.spectra.db` exists, else an indexed `MzMLReader`.
- Spectrum cache is fail-closed: `info`/`is_stale` are source-blind, `verify` refuses a
  network-share source without an explicit local override and exits non-zero on mismatch.
- Continuous integration (GitHub Actions, offline pytest on Python 3.10 and 3.13).
- `CITATION.cff`.

### Changed
- Minimum supported Python is now 3.10. Python 3.9 is end-of-life, and pyteomics ≥ 5.0 (the
  mzML backend) no longer imports on it — caught by the first CI run.
- The structure renderer refuses to guess a topology from a glycan composition.
- Main-only git workflow: pre-commit branch guard and `scripts/dev-setup.sh`.

### Fixed
- XIC area under numpy ≥ 2 (the `trapz` fallback was evaluated eagerly).
- `pyteomics[xml]` is declared as the dependency; a plain `pyteomics` install could not import the package.
- README development-install clone URL pointed at a non-existent GitHub user (#5).

## [0.1.0] - 2026-02-11

### Added
- Initial release: indexed mzML reading, ion search, theoretical fragment calculation and
  peak matching, HCD/EThcD scan pairing, deisotoping, spectral similarity and averaging,
  extracted ion chromatograms, protease definitions.

[Unreleased]: https://github.com/lfu46/mzml-utils/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/lfu46/mzml-utils/releases/tag/v0.1.0
