#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/cq/generated"
w="$PWD/ce-output/pristine"
cb="$w/cr_int_t1d_cb"
cp custody/cq/cq_compare.py "$w/"
python custody/cq/trial_prepare.py | tee "$out/trial-prepare.log"
sha256sum "$cb/"{odd_merged_corrector.npz,odd_combined_projection.npz,odd_merged_tail.npz,odd_BQ_energy.npz} > "$out/merged-inputs.sha256"
(cd "$cb" && python bw_propose_full.py) | tee "$out/proposal.log"
cp "$cb/bw_full_proposal.npz" "$out/proposal16.npz"
for dim in 14 15 16; do
 trial="$out/trial$dim"; mkdir -p "$trial"
 python - "$out/proposal16.npz" "$cb/bw_full_proposal.npz" "$dim" <<'PY'
import sys,numpy as np
a=np.load(sys.argv[1]);n=int(sys.argv[3]);d={k:a[k] for k in a.files}
assert d['Q'].shape==(119,16) and d['X'].shape==(112,119)
for k in ['Q','v0','g0','t','T']:
    assert d[k].ndim==2 and d[k].shape[1]==16,k
    d[k]=d[k][:,:n]
np.savez_compressed(sys.argv[2],**d)
PY
 cp "$cb/bw_full_proposal.npz" "$trial/"
 python "$cb/bw_prepare.py" | tee "$trial/prepare.log"
 cp "$cb/"{bs_trial.npz,odd_joint_proposal.npz} "$trial/"
 for part in 0 1 2 3; do
  python "$cb/bs_circle.py" 32768 "$part" | tee "$trial/circle$part.log"
  cp "$cb/circle_part$part.npz" "$trial/"
 done
 python "$cb/bs_combine.py" 32768 | tee "$trial/combine.log"
 cp "$cb/bs_circle.npz" "$trial/"
 python "$cb/bv_joint_action.py" odd 4096 | tee "$trial/action.log"
 cp "$cb/odd_joint_Uaction.npz" "$trial/"
 python "$w/cq_compare.py" odd "$trial" | tee "$trial/compare.log"
 sha256sum -c "$out/merged-inputs.sha256"
done
python - "$out" <<'PY'
import sys,json,numpy as np
from pathlib import Path
from fractions import Fraction
p=Path(sys.argv[1]);v=np.load(p/'proposal16.npz');ev=np.linalg.eigvalsh(v['central'])
res={}
for n in [14,15,16]:
 a=json.loads((p/f'trial{n}/cj-summary.json').read_text());b=max(a,key=lambda x:Fraction(x['lower']))
 res[str(n)]=dict(physical_correctors=112,dimension=n,proposal_eigenvalue=float(ev[-n]),best=b)
(p/'result.json').write_text(json.dumps(res,indent=2)+'\n')
print('CQ RESULT',json.dumps({k:dict(proposal=v['proposal_eigenvalue'],lower=v['best']['lower'],positive=v['best']['positive']) for k,v in res.items()}),flush=True)
PY
