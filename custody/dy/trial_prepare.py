from pathlib import Path
import shutil,json
import numpy as np
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';cb=w/'cr_int_t1d_cb'
# Authenticate the BR production cross against the untouched CB inputs,
# then apply the known CK proposal-only changes.
src=(r/'custody/ck/prepare.py').read_text().replace("custody/ck/generated","custody/dy/generated")
exec(compile(src,str(r/'custody/ck/prepare.py'),'exec'),dict(__file__=str(r/'custody/ck/prepare.py')))
# Attach the frozen DV physical bank, without changing its inputs.
names=['odd_merged_corrector.npz','odd_combined_projection.npz','odd_merged_tail.npz','odd_BQ_energy.npz']
for n in names:shutil.copyfile(r/'custody/dv/generated/merged'/n,cb/n)
a=np.load(cb/'odd_merged_corrector.npz');b=np.load(cb/'odd_combined_projection.npz');e=np.load(cb/'odd_BQ_energy.npz')
assert a['coefficient'].shape==(491,352)
assert b['Y'].shape==(8192,471) and int(b['d0'])==119
assert e['K'].shape==(352,119) and e['D'].shape==(352,352)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('d=119;m=64')==1 and s.count('ntrial=13')==1 and s.count('fac=1.0;')==1
s=s.replace('d=119;m=64','d=119;m=352').replace('ntrial=13','ntrial=31');s=s.replace('np.linalg.eigvalsh(M)[-16:]','np.linalg.eigvalsh(M)[-36:]');assert s.count('[-4096:]')==1
s=s.replace('[-4096:]','[-8192:]');p.write_text(s)
shutil.copyfile(p,r/'custody/dy/generated/bw_propose_full.py')
print('DV 352 CORRECTOR INPUT SHAPES PASS',flush=True)

# Reuse frozen trial arrays, not a new numerical proposal.
import hashlib
frozen=r/'custody/dv/generated/trial31'
names=['bw_full_proposal.npz','bs_trial.npz','odd_joint_proposal.npz','bs_circle.npz','odd_joint_Uaction.npz']
records={}
for name in names:
    src=frozen/name
    dst=cb/name
    shutil.copyfile(src,dst)
    h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert h(src)==h(dst)
    records[name]=h(dst)
assert np.load(cb/'bs_trial.npz')['Q'].shape==(119,31)
(r/'custody/dy/generated/frozen-trial.json').write_text(json.dumps(records,indent=2)+'\n')
print('DY FROZEN DV TRIAL ATTACHED',flush=True)
