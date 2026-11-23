#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/audits/canonical_return_t1b_exact_audit.py
python3 src/audits/canonical_return_t1c_rational_enclosures.py
python3 src/audits/canonical_return_t1d_uniform_profile.py
