#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# One-time setup: link the private DATA and heavy MODEL artifacts into this repo.
#
# This repo holds CODE only. Two things live OUTSIDE git, in an Insync-synced
# "store" folder, and are symlinked back in by this script:
#   • data/   — the raw whale/gull records (private; the authors do not share them)
#   • models/ — heavy generated objects: posterior sample arrays, prediction
#               samples, z-samples, fitted probabilities, attack p_sim, plus any
#               object that embeds the raw observations (stan data, data subsets,
#               imputed-behaviour files, …). ~7 GB, not for git.
#
# The store location differs per machine, so it is saved to a gitignored file
# (.local-paths) instead of being hardcoded anywhere committed.
#
# USAGE
#   ./setup.sh /path/to/srw_behaviour-store   # first time: give the store path
#   ./setup.sh                                # later: reuses the saved path
#
# See README for where to obtain the store (request the data from the authors).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"                       # always operate from the repo root

# 1. Resolve the machine-local store path. Load anything saved before
#    (gitignored), then let a command-line argument override it.
[ -f .local-paths ] && source .local-paths
[ "$#" -ge 1 ] && STORE_ROOT="$1"

if [ -z "${STORE_ROOT:-}" ]; then
  echo "Usage: ./setup.sh /path/to/srw_behaviour-store"
  echo "(the store holds the private data + heavy model objects — see README)"
  exit 1
fi

# 2. Make sure the store actually exists (the #1 mistake is a wrong path).
if [ ! -d "$STORE_ROOT" ]; then
  echo "ERROR: store folder not found: $STORE_ROOT"
  echo "Obtain it first (see README), then pass the correct path."
  exit 1
fi

# 3. Persist the path for next time, so a later re-run is just "./setup.sh".
echo "STORE_ROOT=$STORE_ROOT" > .local-paths

# 4. Link each heavy folder. The store mirrors the repo's paths (data/, models/),
#    so this is a mechanical loop. Refuse to clobber a real directory.
for rel in data models; do
  target="$STORE_ROOT/$rel"
  if [ -e "$rel" ] && [ ! -L "$rel" ]; then
    echo "ERROR: $rel exists and is not a symlink — refusing to overwrite"; exit 1
  fi
  [ -d "$target" ] || { echo "  note: $rel is empty in the store — creating it"; mkdir -p "$target"; }
  ln -sfn "$target" "$rel"
  echo "linked  $rel  ->  $target"
done

echo
echo "Setup complete — private data and heavy objects are linked into the repo."
echo "Check it worked:  ls data/   &&   ls models/"
