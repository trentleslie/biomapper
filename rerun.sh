#!/usr/bin/env bash
set -euo pipefail
exec /home/trentleslie/agy-team-competition/harmonize/bin/bm-python \
  /home/trentleslie/agy-team-competition/harmonize/winner/scripts/run_dataset_b.py \
  "$@"
