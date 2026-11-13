#!/usr/bin/env bash
# Reproduce all figures and results of the paper (Linux/macOS).
# Usage: ./run_all.sh          full run
#        ./run_all.sh --quick  3 trials per point, output in quick_check/
set -e
cd "$(dirname "$0")"
python3 reproduce.py "$@"
