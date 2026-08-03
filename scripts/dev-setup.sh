#!/bin/sh
# One-time setup after cloning. Activates the versioned pre-commit gate.
#
# core.hooksPath is per-clone local config, not something git can carry in the repo, so a fresh
# clone has no gate until this runs. That is the one manual step.

set -e
REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

git config core.hooksPath .githooks
chmod +x .githooks/pre-commit scripts/branch_health.sh

echo "✓ pre-commit gate activated  (git config core.hooksPath = .githooks)"
echo "  It blocks a commit while any branch other than main exists."
echo
echo "Next: pip install -e .    # mzml_utils importable from any directory"
