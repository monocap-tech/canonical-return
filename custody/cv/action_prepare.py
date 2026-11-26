from pathlib import Path
import zipfile,shutil,hashlib
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';old=w/'cr_int_t1d_ca_odd';new=w/'cr_int_t1d_cv_new';new.mkdir(exist_ok=False)
bo=w/'cr_int_t1d_bo';bo.mkdir(exist_ok=True)
with zipfile.ZipFile(r/'ce-output/archives/CR-INT-T1D-BO-reproducibility.zip') as z:
    for name in ['odd_joint.npz','odd_overlap.npz']:
        member='cr_int_t1d_bo/'+name
        (bo/name).write_bytes(z.read(member))
        print('BO source',member,hashlib.sha256((bo/name).read_bytes()).hexdigest(),flush=True)
for n in ['bn_energy.py','bq_action.py','combine.py','by_action_check.py','by_energy.py']:
    shutil.copyfile(old/n,new/n)
for n in ['odd_corrector.npz','odd_U_tail.npz','odd_new.npy']:
    shutil.copyfile(r/'custody/cv/generated'/n,new/n)
p=new/'by_action_check.py';s=p.read_text()
for a,b in [('F(1,10**54)','F(1,100)'),('F(4,10**55)','F(1,10**30)'),('F(132,10**9)','F(1,10**5)*4')]:
    assert s.count(a)==1
    s=s.replace(a,b)
p.write_text(s)
print('CV source preparation PASS',flush=True)
