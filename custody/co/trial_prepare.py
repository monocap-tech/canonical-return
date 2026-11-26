from pathlib import Path
import shutil,json
import numpy as np
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';cb=w/'cr_int_t1d_cb'
# Authenticate the BR production cross against the untouched CB inputs,
# then apply the known CK proposal-only changes.
src=(r/'custody/ck/prepare.py').read_text().replace("custody/ck/generated","custody/co/generated")
exec(compile(src,str(r/'custody/ck/prepare.py'),'exec'),dict(__file__=str(r/'custody/ck/prepare.py')))
# Attach the independently audited CM merged physical bank.
names=['odd_merged_corrector.npz','odd_combined_projection.npz','odd_merged_tail.npz','odd_BQ_energy.npz']
for n in names:shutil.copyfile(r/'custody/co/generated/merged'/n,cb/n)
a=np.load(cb/'odd_merged_corrector.npz');b=np.load(cb/'odd_combined_projection.npz');e=np.load(cb/'odd_BQ_energy.npz')
assert a['coefficient'].shape==(235,96)
assert b['Y'].shape==(8192,215) and int(b['d0'])==119
assert e['K'].shape==(96,119) and e['D'].shape==(96,96)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('d=119;m=64')==1 and s.count('ntrial=13')==1 and s.count('fac=1.0;')==1
s=s.replace('d=119;m=64','d=119;m=96');p.write_text(s)
shutil.copyfile(p,r/'custody/co/generated/bw_propose_full.py')
print('CO 96 CORRECTOR INPUT SHAPES PASS',flush=True)
