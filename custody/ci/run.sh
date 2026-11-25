#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ci/generated"
w="$PWD/ce-output/pristine"
ca="$w/cr_int_t1d_ca"
frozen="$PWD/custody/ch/generated/expanded-8192-8192/cf-baseline.json"
mkdir -p "$out"
cp custody/ci/ci_compare.py "$w/"
cp custody/cf/generated/{bu_trial.npz,joint_circle.npz} "$ca/"
cp custody/ch/generated/expanded-8192-8192/{even_joint_proposal.npz,even_joint_Uaction.npz} "$ca/"
sha256sum "$ca/bu_trial.npz" "$ca/joint_circle.npz" > "$out/frozen-physical.sha256"
compare() {
 local tag="$1"; shift
 mkdir -p "$out/$tag"
 python "$w/ci_compare.py" even "$out/$tag" "$frozen" "$@" | tee "$out/$tag/compare.log"
 sha256sum -c "$out/frozen-physical.sha256"
}
compare control-8192-8192 --control
trial() {
 local tag="$1" exterior="$2" interior="$3"
 mkdir -p "$out/$tag"
 python - "$ca/even_joint_proposal.npz" "$interior" <<'PY'
import numpy as np,sys
a=np.load(sys.argv[1]);assert len(a['j'])==int(sys.argv[2]) and a['T'].shape==(int(sys.argv[2]),4)
PY
 cp "$ca/even_joint_proposal.npz" "$out/$tag/"
 python "$ca/bv_joint_action.py" even "$exterior" | tee "$out/$tag/action.log"
 cp "$ca/even_joint_Uaction.npz" "$out/$tag/"
 compare "$tag"
}
trial matched-8192-16384 16384 8192
python custody/ci/propose.py "$ca" | tee "$out/propose.log"
trial expanded-16384-16384 16384 16384
python - "$out" <<'PY'
import json,sys
from pathlib import Path
from fractions import Fraction
p=Path(sys.argv[1]);res={}
for f in sorted(p.glob('*/summary.json')):
 r=json.loads(f.read_text());best=max(r.values(),key=lambda v:Fraction(v['lower']));res[f.parent.name]=best
(p/'results.json').write_text(json.dumps(res,indent=2)+'\n')
print('CI BEST RESULTS',json.dumps(res),flush=True)
PY
