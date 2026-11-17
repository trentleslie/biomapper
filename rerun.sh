#!/usr/bin/env bash
set -euo pipefail
exec /home/trentleslie/agy-team-competition/harmonize/bin/bm-python \
  --venv /home/trentleslie/agy-team-competition/harmonize/worktrees/team-a/.venv \
  /home/trentleslie/agy-team-competition/harmonize/worktrees/team-a/scripts/run_dataset_b.py \
  "$@"
