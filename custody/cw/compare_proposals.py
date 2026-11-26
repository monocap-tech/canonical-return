from pathlib import Path
import numpy as np,json
r=Path(__file__).resolve().parents[2]
a=np.load(r/'custody/cv/generated/proposal17.npz');b=np.load(r/'custody/cw/generated/proposal17.npz')
assert a['j'].shape==(4096,) and b['j'].shape==(8192,)
assert len(np.unique(b['j']))==8192 and np.all(b['j']%2==1)
assert np.array_equal(np.intersect1d(a['j'],b['j']),a['j'])
keys=['Q','X','v0','g0','t','alpha','central']
match={k:bool(np.array_equal(a[k],b[k])) for k in keys}
# Exact equality is a diagnostic, not a premise of the new directed comparison.
d=dict(old_inner=4096,new_inner=8192,shared_indices=4096,matching_fields=match,physical_bank='CV frozen 160 columns')
print('CW PROPOSAL COMPARISON',json.dumps(d),flush=True)
(r/'custody/cw/generated/proposal-comparison.json').write_text(json.dumps(d,indent=2)+'\n')
