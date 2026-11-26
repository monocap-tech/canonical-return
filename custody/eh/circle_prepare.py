"""Prepare a matched circular replay at the old signed-symbol cutoff.
The full physical representative includes high Neumann modes, so all source
labels must be included in the Gram evaluation even outside the circle cutoff.
"""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
root=Path(__file__).resolve().parents[2]
w=root/'ce-output/pristine';cb=w/'cr_int_t1d_cb';out=root/'custody/eh/generated';dz=root/'custody/dz/generated/trial31'
a=np.load(cb/'bs_trial.npz')
assert a['yn'].shape==(32768,31)
C=32768
jall=np.r_[a['js'],2*a['ms']]
assert len(set(jall))==len(jall) and min(jall)>0
ids=np.union1d(jall,np.arange(1,C,2,dtype=np.int64))
parts=(len(ids)+8191)//8192
symbol=np.load(w/'cr_int_t1d_bs/bs_symbol.npz')
assert np.array_equal(symbol['ids'],np.arange(1,C,2))
# Same inverse response and circular symbol rows; original inverse-action
# producer does not depend on the finite representative except for these rows.
u=np.load(dz/'odd_joint_Uaction.npz');p=np.load(dz/'odd_joint_proposal.npz')
assert np.array_equal(u['ids'],symbol['ids']) and np.array_equal(u['j'],p['j']) and np.array_equal(u['T'],p['T'])
assert np.array_equal(np.load(cb/'odd_joint_proposal.npz')['T'],p['T'])
for name in ('bs_circle.py','bs_combine.py'):
    src=root/'custody/ee/generated/sources/cr_int_t1d_cb'/name
    dest=cb/name;assert src.read_bytes()==dest.read_bytes()
    s=src.read_text()
    # Reduce row tile width for memory only, not the logical 8192-row part.
    for old,new in [('range(startrow,endrow,256)','range(startrow,endrow,32)'),('min(b+256,endrow)','min(b+32,endrow)')]:
        assert s.count(old)==1;s=s.replace(old,new)
    dest.write_text(s);shutil.copyfile(dest,out/name)
# exact matched inverse action custody: no new inverse coefficient/row.
shutil.copyfile(dz/'odd_joint_Uaction.npz',cb/'odd_joint_Uaction.npz')
r=dict(physical_cutoff=65536,circular_cutoff=C,source_labels=len(jall),max_source_label=int(max(jall)),all_gram_rows=len(ids),parts=parts,inverse_action_sha256=hashlib.sha256((cb/'odd_joint_Uaction.npz').read_bytes()).hexdigest(),scope='matched circle cutoff; full finite-representative Gram required; omitted circular tail still enclosed')
(out/'circle-plan.json').write_text(json.dumps(r,indent=2)+'\n')
print(parts)
