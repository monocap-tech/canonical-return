#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/cf/generated"
mkdir -p "$out"
w="$PWD/ce-output/trial4"
ca="$w/cr_int_t1d_ca"
python custody/ce/prepare_four.py | tee "$out/prepare-scratch.log"
python "$ca/bx_prepare.py" | tee "$out/prepare.log"
for n in 0 1 2; do python "$ca/bu_circle.py" "$n" | tee "$out/circle-$n.log"; done
python "$ca/bx_combine.py" | tee "$out/combine.log"
cp "$ca"/{bx_prepare.py,bu_trial.npz,joint_circle.npz,raw_part0.npz,raw_part1.npz,raw_part2.npz} "$out/"
cp "$w/ce_four.py" "$out/"
# Freeze the physical representation throughout all matched trials.
sha256sum "$ca/bu_trial.npz" "$ca/joint_circle.npz" > "$out/frozen-physical.sha256"
check() {
 local tag="$1" exterior="$2"
 mkdir -p "$out/$tag"
 cp "$ca/even_joint_proposal.npz" "$out/$tag/"
 python "$ca/bv_joint_action.py" even "$exterior" | tee "$out/$tag/action.log"
 cp "$ca/even_joint_Uaction.npz" "$out/$tag/"
 python "$w/ce_four.py" even | tee "$out/$tag/comparison.log"
 cp "$w/cd_results/even.json" "$out/$tag/comparison.json"
 sha256sum -c "$out/frozen-physical.sha256"
}
check original-2048-2048 2048
python - "$out/original-2048-2048/comparison.json" <<'PY'
import sys,json
from fractions import Fraction
r=json.load(open(sys.argv[1]));assert r['n']==4 and Fraction(r['lower'])==Fraction(-7,1250)
print('CE BASELINE EXACT MATCH',flush=True)
PY
check matched-2048-4096 4096
# New numerical proposal only; same exact physical v, t, q and x.
python custody/cf/propose.py "$ca" | tee "$out/propose.log"
check expanded-4096-4096 4096
python - "$out" <<'PY'
import sys,json
from pathlib import Path
p=Path(sys.argv[1]);results={}
for q in sorted(p.glob('*/comparison.json')):
 r=json.loads(q.read_text());results[q.parent.name]=r
print('CF EXACT RESULTS',json.dumps(results),flush=True)
(p/'results.json').write_text(json.dumps(results,indent=2)+'\n')
PY
