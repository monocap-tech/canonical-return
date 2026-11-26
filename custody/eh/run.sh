#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/eh/generated"
cb="$PWD/ce-output/pristine/cr_int_t1d_cb"
python custody/eh/synthesis.py | tee "$out/synthesis.log"
parts="$(python custody/eh/circle_prepare.py)"
for ((part=0;part<parts;part++)); do
 python "$cb/bs_circle.py" 32768 "$part" | tee "$out/circle-$part.log"
 cp "$cb/circle_part$part.npz" "$out/"
done
python "$cb/bs_combine.py" 32768 | tee "$out/combine.log"
cp "$cb/bs_circle.npz" "$out/"
python - "$out" <<'PY'
import sys,json,numpy as np
from pathlib import Path
p=Path(sys.argv[1]);a=np.load(p/'bs_circle.npz')
assert np.array_equal(a['ids'],np.arange(1,32768,2))
assert a['Z'].shape==(16384,31) and a['G'].shape==(31,31)
assert len(a['reprerr'])==31
print('EH CIRCLE COMPLETE',a['Z'].shape,flush=True)
PY
