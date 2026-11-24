#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 - <<'PY'
from pathlib import Path
import hashlib
p=Path('CR-ALP-034_initial_rational_intervals_1e25.txt')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2'
print('PASS: historical source input SHA256 custody')
PY
g++ -O2 -std=c++17 cr_alp_034_midpoint_12_mode_interval_core.cpp -o cr_alp_034_midpoint_core
./cr_alp_034_midpoint_core CR-ALP-034_initial_rational_intervals_1e25.txt 30000 > CR-ALP-034_midpoint_core.log
python3 cr_alp_034_independent_midpoint_ldl_verifier.py
