#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 cr_alp_038_combined_infinite_tail_rational_audit.py
python3 cr_alp_038_full_interior_infinite_schur_audit.py
