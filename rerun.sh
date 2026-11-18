#!/usr/bin/env bash
set -euo pipefail
exec "/home/trentleslie/agy-team-competition/harmonize/bin/bm-python" --venv "/home/trentleslie/agy-team-competition/harmonize/worktrees/team-c/.venv" "/home/trentleslie/agy-team-competition/harmonize/worktrees/team-c/scripts/run_team_c.py" "$@"
