#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 - <<'PY'
from pathlib import Path
from hashlib import sha256
checks={
'CR-ALP-040_initial_rational_intervals_1e25.txt':'6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9',
'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv':'7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c',
}
for file,hash in checks.items():
 assert sha256(Path(file).read_bytes()).hexdigest()==hash
 print('PASS pinned source interval SHA256',file)
PY
python3 cr_alp_040_shell_geometry_audit.py
python3 cr_alp_040_stable_14_mode_sylvester_audit.py
python3 cr_alp_040_mixed_odd_even_shell_audit.py
printf '\nCR-ALP-040 FAST EXACT SOURCE BANK REPRODUCTION PASS\n'
