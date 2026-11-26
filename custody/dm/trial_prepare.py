from pathlib import Path
import shutil,json
import numpy as np
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';cb=w/'cr_int_t1d_cb'
# Authenticate the BR production cross against the untouched CB inputs,
# then apply the known CK proposal-only changes.
src=(r/'custody/ck/prepare.py').read_text().replace("custody/ck/generated","custody/dm/generated")
exec(compile(src,str(r/'custody/ck/prepare.py'),'exec'),dict(__file__=str(r/'custody/ck/prepare.py')))
# Attach the frozen DK physical bank, without changing its inputs.
names=['odd_merged_corrector.npz','odd_combined_projection.npz','odd_merged_tail.npz','odd_BQ_energy.npz']
for n in names:shutil.copyfile(r/'custody/dk/generated/merged'/n,cb/n)
a=np.load(cb/'odd_merged_corrector.npz');b=np.load(cb/'odd_combined_projection.npz');e=np.load(cb/'odd_BQ_energy.npz')
assert a['coefficient'].shape==(427,288)
assert b['Y'].shape==(8192,407) and int(b['d0'])==119
assert e['K'].shape==(288,119) and e['D'].shape==(288,288)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('d=119;m=64')==1 and s.count('ntrial=13')==1 and s.count('fac=1.0;')==1
s=s.replace('d=119;m=64','d=119;m=288').replace('ntrial=13','ntrial=26');s=s.replace('np.linalg.eigvalsh(M)[-16:]','np.linalg.eigvalsh(M)[-26:]');assert s.count('[-4096:]')==1
s=s.replace('[-4096:]','[-8192:]');p.write_text(s)
shutil.copyfile(p,r/'custody/dm/generated/bw_propose_full.py')
print('DK 288 CORRECTOR INPUT SHAPES PASS',flush=True)

# Reuse frozen trial arrays, not a new numerical proposal.
import hashlib
frozen=r/'custody/dk/generated/trial26'
names=['bw_full_proposal.npz','bs_trial.npz','odd_joint_proposal.npz','bs_circle.npz','odd_joint_Uaction.npz']
records={}
for name in names:
    src=frozen/name
    dst=cb/name
    shutil.copyfile(src,dst)
    h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert h(src)==h(dst)
    records[name]=h(dst)
assert np.load(cb/'bs_trial.npz')['Q'].shape==(119,26)
(r/'custody/dm/generated/frozen-trial.json').write_text(json.dumps(records,indent=2)+'\n')
print('DM FROZEN DK TRIAL ATTACHED',flush=True)
