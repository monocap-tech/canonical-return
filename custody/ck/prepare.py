from pathlib import Path
import zipfile,shutil,hashlib
import numpy as np
r=Path(__file__).resolve().parents[2]; w=r/'ce-output/pristine'; cb=w/'cr_int_t1d_cb'; ca=w/'cr_int_t1d_ca_odd'
# Use the literal production path, not the pilot4096 checkpoint.
import io,json
archive=r/'ce-output/archives/CR-INT-T1D-BR-reproducibility.zip'
with zipfile.ZipFile(archive) as z:
    root=z.read('odd_surplus_cross.npz')
    pilot=z.read('pilot4096/odd_surplus_cross.npz')
old=np.load(cb/'bw_full_proposal.npz')
trial=np.load(cb/'bs_trial.npz')
prod=np.load(io.BytesIO(root))
assert np.array_equal(prod['weights'],old['alpha'])
assert np.array_equal(prod['weights'],trial['a'])
assert len(prod['indices'])==old['t'].shape[0]==1882
assert not np.intersect1d(prod['indices'],np.load(w/'cr_int_t1d_bk/odd/physical_dictionary.npz')['indices']).size
br=w/'cr_int_t1d_br';br.mkdir(exist_ok=True)
(br/'odd_surplus_cross.npz').write_bytes(root)
out=r/'custody/ck/generated';out.mkdir(exist_ok=True)
record=dict(archive=archive.name,member='odd_surplus_cross.npz',sha256=hashlib.sha256(root).hexdigest(),pilot_sha256=hashlib.sha256(pilot).hexdigest(),production_weights_match_both=True,rows=len(prod['indices']))
(out/'br-provenance.json').write_text(json.dumps(record,indent=2)+'\n')
print('BR PRODUCTION CROSS',record,flush=True)
for n in ['bw_prepare.py','bs_circle.py','bs_combine.py']:
    shutil.copyfile(ca/n,cb/n)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('fac=1.02;')==1 and s.count('ntrial=min(16,max(9,int(np.sum(ev>0))))')==1
s=s.replace('fac=1.02;','fac=1.0;').replace('ntrial=min(16,max(9,int(np.sum(ev>0))))','ntrial=13')
p.write_text(s)
print('THIRTEEN COLUMN PROPOSAL SOURCE READY',flush=True)
