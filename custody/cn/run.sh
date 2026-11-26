#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/cn/generated"
w="$PWD/ce-output/pristine"
cb="$w/cr_int_t1d_cb"
mkdir -p "$out/trial13"
cp custody/cn/cn_compare.py "$w/"
python custody/cn/prepare.py | tee "$out/prepare.log"
sha256sum "$cb/"{odd_merged_corrector.npz,odd_combined_projection.npz,odd_merged_tail.npz,odd_BQ_energy.npz} > "$out/merged-inputs.sha256"
(cd "$cb" && python bw_propose_full.py) | tee "$out/trial13/proposal.log"
cp "$cb/bw_full_proposal.npz" "$out/trial13/"
python "$cb/bw_prepare.py" | tee "$out/trial13/prepare.log"
cp "$cb/"{bs_trial.npz,odd_joint_proposal.npz} "$out/trial13/"
python "$cb/bs_circle.py" 32768 0 | tee "$out/trial13/circle0.log"
cp "$cb/circle_part0.npz" "$out/trial13/"
python "$cb/bs_circle.py" 32768 1 | tee "$out/trial13/circle1.log"
cp "$cb/circle_part1.npz" "$out/trial13/"
python "$cb/bs_circle.py" 32768 2 | tee "$out/trial13/circle2.log"
cp "$cb/circle_part2.npz" "$out/trial13/"
python "$cb/bs_circle.py" 32768 3 | tee "$out/trial13/circle3.log"
cp "$cb/circle_part3.npz" "$out/trial13/"
python "$cb/bs_combine.py" 32768 | tee "$out/trial13/combine.log"
cp "$cb/bs_circle.npz" "$out/trial13/"
python "$cb/bv_joint_action.py" odd 4096 | tee "$out/trial13/action.log"
cp "$cb/odd_joint_Uaction.npz" "$out/trial13/"
python "$w/cn_compare.py" odd "$out/trial13" | tee "$out/trial13/compare.log"
cp "$cb/"{bw_full_proposal.npz,bs_trial.npz,bs_circle.npz,odd_joint_proposal.npz,odd_joint_Uaction.npz,circle_part0.npz,circle_part1.npz,circle_part2.npz,circle_part3.npz} "$out/trial13/"

sha256sum -c "$out/merged-inputs.sha256"
python - "$out" <<'PY'
import sys,json,numpy as np
from pathlib import Path
from fractions import Fraction
p=Path(sys.argv[1]);v=np.load(p/'trial13/bw_full_proposal.npz')
assert v['Q'].shape==(119,13) and v['X'].shape==(80,119)
a=json.loads((p/'trial13/cj-summary.json').read_text())
best=max(a,key=lambda x:Fraction(x['lower']))
d=dict(physical_correctors=80,dimension=13,proposal_eigenvalue=float(np.linalg.eigvalsh(v['central'])[-13]),best=best)
(p/'result.json').write_text(json.dumps(d,indent=2)+'\n')
print('CN RESULT',json.dumps(d),flush=True)
PY
