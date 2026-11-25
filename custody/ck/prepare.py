from pathlib import Path
import zipfile,shutil,hashlib
import numpy as np
r=Path(__file__).resolve().parents[2]; w=r/'ce-output/pristine'; cb=w/'cr_int_t1d_cb'; ca=w/'cr_int_t1d_ca_odd'
# Locate BR's exact archived cross checkpoint, requiring identical bytes for duplicates.
hits=[]
for p in sorted((r/'ce-output/archives').glob('*.zip')):
    with zipfile.ZipFile(p) as z:
        for n in z.namelist():
            if Path(n).name=='odd_surplus_cross.npz':
                hits.append((p.name,n,z.read(n)))
assert hits, 'missing BR odd cross'
assert len({hashlib.sha256(x[2]).hexdigest() for x in hits})==1,[(x[0],x[1]) for x in hits]
br=w/'cr_int_t1d_br';br.mkdir(exist_ok=True)
(br/'odd_surplus_cross.npz').write_bytes(hits[0][2])
print('BR cross provenance',[(x[0],x[1]) for x in hits],hashlib.sha256(hits[0][2]).hexdigest(),flush=True)
for n in ['bw_prepare.py','bs_circle.py','bs_combine.py']:
    shutil.copyfile(ca/n,cb/n)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('fac=1.02;')==1 and s.count('ntrial=min(16,max(9,int(np.sum(ev>0))))')==1
s=s.replace('fac=1.02;','fac=1.0;').replace('ntrial=min(16,max(9,int(np.sum(ev>0))))','ntrial=13')
p.write_text(s)
print('THIRTEEN COLUMN PROPOSAL SOURCE READY',flush=True)
