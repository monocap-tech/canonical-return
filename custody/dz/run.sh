#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/dz/generated"
w="$PWD/ce-output/pristine"
cb="$w/cr_int_t1d_cb"
trial="$out/trial31"
mkdir -p "$trial" "$out/block"
cp custody/dz/dz_compare.py custody/dz/block_compare.py "$w/"
python custody/dz/trial_prepare.py | tee "$out/trial-prepare.log"
sha256sum "$cb/"{odd_merged_corrector.npz,odd_combined_projection.npz,odd_merged_tail.npz,odd_BQ_energy.npz} > "$out/merged-inputs.sha256"
python custody/dz/propose.py | tee "$out/proposal.log"
cp "$cb/bw_full_proposal.npz" "$trial/"
python "$cb/bw_prepare.py" | tee "$trial/prepare.log"
cp "$cb/"{bs_trial.npz,odd_joint_proposal.npz} "$trial/"
for part in 0 1 2 3; do
 python "$cb/bs_circle.py" 32768 "$part" | tee "$trial/circle$part.log"
 cp "$cb/circle_part$part.npz" "$trial/"
done
python "$cb/bs_combine.py" 32768 | tee "$trial/combine.log"
cp "$cb/bs_circle.npz" "$trial/"
python "$cb/bv_joint_action.py" odd 8192 | tee "$trial/action.log"
cp "$cb/odd_joint_Uaction.npz" "$trial/"
python "$w/dz_compare.py" odd "$trial" | tee "$trial/compare.log"
python "$w/block_compare.py" odd "$out/block" "$trial" | tee "$out/block/compare.log"
sha256sum -c "$out/merged-inputs.sha256"
python - "$out" <<'PY'
import json,sys
from pathlib import Path
from fractions import Fraction as F
p=Path(sys.argv[1]);a=json.loads((p/'trial31/cj-summary.json').read_text());b=json.loads((p/'block/result.json').read_text());scalar=max(a,key=lambda x:F(x['lower']));best=max([scalar,b['best']],key=lambda x:F(x['lower']))
r=dict(dimension=31,physical_correctors=352,proposal=json.loads((p/'proposal-check.json').read_text()),scalar_best=scalar,block_best=b['best'],best=best)
(p/'result.json').write_text(json.dumps(r,indent=2)+'\n')
print('DZ BEST',best['lower'],best['positive'],flush=True)
PY
