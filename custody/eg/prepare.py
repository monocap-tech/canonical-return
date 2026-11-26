from pathlib import Path
import json,shutil,hashlib
import numpy as np
root=Path(__file__).resolve().parents[2];w=root/'ce-output/pristine';out=root/'custody/eg/generated'
p=w/'cr_int_t1d_eg_physical';p.mkdir(exist_ok=False)
h=w/'cr_int_t1d_eg_head';h.mkdir(exist_ok=False)
src=root/'custody/ee/generated/sources/cr_int_t1d_ca_odd/bq_action.py'
assert src.read_bytes()==(w/'cr_int_t1d_ca_odd/bq_action.py').read_bytes()
shutil.copyfile(src,p/'bq_action.py')
for folder in [p,h]:shutil.copyfile(w/'cr_int_t1d_ca_odd/bn_energy.py',folder/'bn_energy.py')
shutil.copyfile(root/'custody/dv/generated/merged/odd_merged_corrector.npz',p/'odd_corrector.npz')
assert np.load(p/'odd_corrector.npz')['coefficient'].shape==(491,352)
source=root/'custody/ef/generated/head-sources/BP/bp_head.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='bfe0ea71b464b2ea6393579e739e5f4a42670c6f1cc1dc62af50ea4a70923a83'
s=source.read_text()
old="h=np.load(bo/(par+'_head.npz'));headindices=list(map(int,h['indices']));assert ms[:len(headindices)]==headindices;old=obj(h['action']);olds=int(h['scale']);Yold=(old@B*YS)//(olds*BS);err_old=up(norm(old,olds)*eB+F(str(h['action_error']))*(nB+eB)+F(Yold.size,YS));print(par,'old normalized head error',float(err_old),flush=True)"
new="h=np.load(root.parent/'cr_int_t1d_bp'/(par+'_head_full.npz'));headindices=list(map(int,h['indices']));assert ms[:len(headindices)]==headindices and len(headindices)==8192 and int(h['scale'])==YS;Yold=h['Y'];err_old=F(str(h['error']));assert Yold.shape==(8192,119);print(par,'frozen BP normalized head error',float(err_old),flush=True)"
assert s.count(old)==1;s=s.replace(old,new)
# The arithmetic formula and directed primitive remain unchanged; only prefix reuse and additional high-frequency EM overlap probes change.
anchor='[258+p,510+p,1022+p,4094+p,8190+p,16382+p]'
assert s.count(anchor)==1;s=s.replace(anchor,'[258+p,510+p,1022+p,4094+p,8190+p,16382+p,32766+p,65534+p]')
(h/'bp_head_extend.py').write_text(s)
(out/'bp_head_extend.py').write_text(s)
# Original BO inputs were verified by reconstruction and recovery.
for name in ['odd_joint.npz','odd_normalization.npz']:
    assert (w/'cr_int_t1d_bo'/name).is_file()
print('EG frozen bank and head source ready',flush=True)
