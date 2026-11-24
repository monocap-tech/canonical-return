#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 cr_alp_040_initial_rational_data.py > CR-ALP-040_initial_generator.log
python3 - <<'PY'
from pathlib import Path
from hashlib import sha256
old = {}
for line in Path('predecessor_source/CR-ALP-034_initial_rational_intervals_1e25.txt').read_text().splitlines():
 k,lo,hi = line.split();old[k] = (lo,hi)
new = {}
for line in Path('CR-ALP-040_initial_rational_intervals_1e25.txt').read_text().splitlines():
 k,lo,hi = line.split();new[k] = (lo,hi)
assert len(old)==177 and len(new)==233
assert all(new.get(k)==v for k,v in old.items())
assert sha256(Path('CR-ALP-040_initial_rational_intervals_1e25.txt').read_bytes()).hexdigest()=='6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9'
print('PASS historical 177 original source-interval entries unchanged; 56 source extension entries added')
PY
g++ -O2 -std=c++17 cr_alp_040_build_signed_14_mode_lower_bank.cpp -o cr_alp_040_bank_builder
./cr_alp_040_bank_builder CR-ALP-040_initial_rational_intervals_1e25.txt 30000 > CR-ALP-040_core_reproduction.log
python3 - <<'PY'
from hashlib import sha256
from pathlib import Path
p=Path('CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv')
assert sha256(p.read_bytes()).hexdigest()=='7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c'
assert len(p.read_text().splitlines())==98
print('PASS 14-mode source-signed midpoint bank SHA-256 and 98 stored entries')
PY
python3 cr_alp_040_shell_geometry_audit.py
python3 cr_alp_040_stable_14_mode_sylvester_audit.py
python3 cr_alp_040_mixed_odd_even_shell_audit.py
printf '\nCR-ALP-040 FULL REPRODUCTION PASS (finite E14 + 13-mode mixed subspace; not RH)\n'
