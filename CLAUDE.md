# CLAUDE.md — mzml-utils

Shared utility library for the Wu lab's MS and glycoproteomics projects. Installed editable
(`pip install -e .`) and imported flat as `mzml_utils` from every other repo.

## Git workflow — `main` only, no branches, no exceptions

**Work on `main`.** `git pull` before starting, `git commit` + `git push` when done.

**Never create a branch or worktree, for a change of any size.** This is enforced, not advised:

- `.githooks/pre-commit` runs `scripts/branch_health.sh` and blocks any commit made while a
  branch other than `main`, or a second worktree, exists. Run `bash scripts/dev-setup.sh` once
  after cloning to activate it (`core.hooksPath` is per-clone and cannot be carried in the repo).
- A global `PreToolUse` hook (`wu-lab-skills/hooks/no_branches.py`) denies the branch-creating
  command outright. Deleting branches and moving between existing refs stay allowed.

**Why this repo in particular.** Every other lab project consumes it through an editable install,
so *what is importable is whatever is checked out on disk*. When `mzml_utils.structure` existed
only on the `p1-cache-safety` branch, a plain `git checkout main` would have silently removed it
from every project on the machine — no error at checkout, just `ModuleNotFoundError` later,
somewhere else. A branch here is not a local decision; it changes what other repos can import.

The rule also has history: NGlycoMM's equivalent guidance once allowed "a branch for a big or
risky change". That judgment call was exercised 16 times and was wrong all 16 — 16 branches, 10
unmerged, 14 of 14 tested pairs conflicting on identical lines. The clause is gone in both repos.

If a branch ever seems genuinely necessary, that is a decision to raise, not a default to take:
`NGLYCO_ALLOW_BRANCH=1` unblocks the hook, and `git commit --no-verify` bypasses the gate.

## Reading spectra

`mzml_utils.open_spectra(path)` is the single entry point — cache-aware, returning a fast local
`SpectrumCache` when one exists beside the mzML, else an indexed `MzMLReader`. Identical
`Spectrum` objects either way. Never use `pyteomics.mzml` directly, and do not instantiate
`MzMLReader` for reading.

## Structure databases

`mzml_utils.structure` is the single place a protein-structure URL is built — never hand-write an
AlphaFold, PDB or UniProt URL. See the `protein-structure` skill for the full contract.

## Tests

```bash
python3 -m pytest tests/ -q -m "not network"   # offline, fast
python3 -m pytest tests/ -q -m network         # hits EBI/UniProt
```

Network-marked tests are anchors against live APIs and are expected to be run deliberately, not
in a tight loop.
