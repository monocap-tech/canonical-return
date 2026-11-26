"""EH: reconstruct the exact finite representative at physical cutoff 65536.
Requires completed, authenticated EG evidence in the current checkout.
No new numerical proposal; does not establish positivity.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
import numpy as np
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parents[2]
eg=root/'custody/eg/generated';dv=root/'custody/dv/generated';dz=root/'custody/dz/generated/trial31'
out=root/'custody/eh/generated';out.mkdir(parents=True,exist_ok=True)
cb=root/'ce-output/pristine/cr_int_t1d_cb'
S=10**12
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def authenticate(p):
    for line in (p/'sha256.txt').read_text().splitlines():
        digest,name=line.split(maxsplit=1)
        assert sha(p/name.removeprefix('./'))==digest,name
for p in (eg,dv,root/'custody/dz/generated'):authenticate(p)
report=json.loads((eg/'extension.json').read_text())
assert report['cutoff']==65536 and report['head_half_target_met'] and report['all_block_half_targets_met']
old=np.load(dz/'bs_trial.npz');prop=np.load(dz/'bw_full_proposal.npz')
for k in ('Q','X','t','v0','g0'):assert np.array_equal(old[k],prop[k]),k
assert old['Q'].shape==(119,31) and old['X'].shape==(352,119)
assert np.array_equal(old['ms'],np.arange(1,16384,2))
assert old['yn'].shape==(8192,31) and int(old['scale'])==S
head=np.load(eg/'odd_head_full.npz');assert head['Y'].shape==(32768,119)
assert np.array_equal(head['indices'],np.arange(1,65536,2)) and int(head['scale'])==S
p=np.load(dv/'merged/odd_combined_projection.npz')
assert np.array_equal(head['Y'][:8192],p['Y'][:,:119])
# Exact same Bq and integer floor as bw_prepare.py. Check the frozen prefix.
Bq=np.vstack([old['Q'].astype(object)*S,old['X'].astype(object)@old['Q'].astype(object)])
prefix=(p['Y'].astype(object)@Bq)//(S*S)
assert np.array_equal(prefix,old['yn'].astype(object))
chunks=[];covered=[];sq=[0]*31
for r in report['new_chunks']:
    lo,hi=r['lo'],r['hi']
    with np.load(eg/f'odd_U_{lo}_{hi}.npz') as a:
        ids=a['indices'];Y=a['Y'];assert np.array_equal(ids,np.arange(lo,hi+1,2))
        assert Y.shape==(len(ids),352) and int(a['scale'])==S
        for k in range(0,len(ids),256):
            ii=ids[k:k+256];Y0=np.column_stack([head['Y'][(ii-1)//2],Y[k:k+256]])
            yq=Y0.astype(object)@Bq;yn=yq//(S*S)
            assert max(abs(v) for v in yn.flat)<2**63
            R=yq-yn*S*S
            for j in range(31):sq[j]+=sum(int(v)**2 for v in R[:,j])
            chunks.append(yn.astype(np.int64));covered.extend(map(int,ii))
assert covered==list(range(16385,65536,2))
yn=np.concatenate([old['yn'],*chunks]);assert yn.shape==(32768,31)
extra=[F(isqrt(v)+(isqrt(v)**2<v),S**3) for v in sq]
errors=[F(str(a))+b for a,b in zip(old['errors'],extra)]
fields={k:old[k] for k in old.files}
fields.update(yn=yn,ms=head['indices'],errors=np.array(list(map(str,errors))))
for k in old.files:
    if k not in ('yn','ms','errors'):assert np.array_equal(fields[k],old[k]),k
np.savez_compressed(cb/'bs_trial.npz',**fields)
np.savez_compressed(out/'bs_trial.npz',**fields)
# New tail/error metadata is separate from the old 8192-row matrix.
# Later receiver must explicitly select it, not silently reinterpret old arrays.
metadata=dict(cutoff=65536,head_error=report['head_error'],head_tail=report['head_tail'],physical_error=report['physical_finite_error'],physical_tail=report['physical_tail'],physical_blocks=report['physical_blocks'])
(out/'representation-budget.json').write_text(json.dumps(metadata,indent=2)+'\n')
summary=dict(scope='finite physical representative only',source_head_sha256=sha(eg/'odd_head_full.npz'),source_frozen_trial_sha256=sha(dz/'bs_trial.npz'),prefix_exact=True,rows=len(yn),columns=31,additional_rows=len(covered),max_new_rounding_allowance=float(max(extra)),max_combined_column_error=float(max(errors)),frozen_fields=[k for k in old.files if k not in ('yn','ms','errors')])
(out/'synthesis.json').write_text(json.dumps(summary,indent=2)+'\n')
print('EH SYNTHESIS',summary,flush=True)
