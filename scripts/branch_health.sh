#!/bin/sh
# branch_health.sh -- block a commit while any branch or extra worktree exists.
#
# Companion to wu-lab-skills/hooks/no_branches.py, which is registered as a global PreToolUse
# hook and denies the branch-creating command outright. That is the better guard because
# nothing gets created at all -- but it only binds Claude Code. It does not bind you in a
# terminal, Codex, an IDE's git UI, or any other tool. This one binds the repository itself.
#
# Why the rule has no size exception. mzml-utils carried `p1-cache-safety` as an unmerged
# feature branch, and an unrelated change (the structure subpackage) landed on top of it. That
# is the shape the NGlycoMM guard was written for after its escape clause -- "branch only for a
# big or risky change" -- was exercised 16 times and was wrong all 16: 16 branches, 10 unmerged,
# 14 of 14 tested pairs conflicting on identical lines.
#
# There is a sharper failure specific to THIS repo: it is consumed by every other lab project
# through `pip install -e`, so what is importable is whatever is checked out on disk. While the
# structure package existed only on a branch, `git checkout main` would have silently removed
# `mzml_utils.structure` from every project on the machine.
#
# Deleting branches is not blocked -- only committing while one exists. Bypass: git commit --no-verify

REPO="$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0
cd "$REPO" || exit 0

rc=0

# --- local branches other than main -------------------------------------------------------
others=$(git for-each-ref --format='%(refname:short)' refs/heads/ | grep -v '^main$' || true)
if [ -n "$others" ]; then
  echo "✗ COMMIT BLOCKED — branches other than main exist:"
  echo "$others" | sed 's/^/    /'
  echo "  No branches on this repo, with no size exception. Land the work on main"
  echo "  (git merge), then: git branch -d <name> && git push origin --delete <name>"
  rc=1
fi

# --- extra worktrees ----------------------------------------------------------------------
n_wt=$(git worktree list --porcelain | grep -c '^worktree ' || true)
if [ "${n_wt:-1}" -gt 1 ]; then
  echo "✗ COMMIT BLOCKED — $n_wt worktrees exist; only the main checkout is allowed:"
  git worktree list | sed 's/^/    /'
  echo "  Remove them: git worktree remove <path>  (then git worktree prune)"
  rc=1
fi

# --- remote branches other than main (warn only) -------------------------------------------
# A warning, not a block: deleting a remote ref needs network and must not gate a local commit.
#
# Filter on the FULL refname, not %(refname:short). git shortens refs/remotes/origin/HEAD to
# plain "origin", which no `origin/...` pattern matches -- the first version of this check
# duly reported "origin" as a stray branch.
remotes=$(git for-each-ref --format='%(refname)' refs/remotes/origin/ 2>/dev/null \
          | grep -v -E '^refs/remotes/origin/(main|HEAD)$' \
          | sed 's|^refs/remotes/||' || true)
if [ -n "$remotes" ]; then
  n=$(echo "$remotes" | wc -l | tr -d ' ')
  echo "  note: $n remote branch(es) still on origin — delete when convenient:"
  echo "$remotes" | sed 's/^/    /'
fi

exit $rc
