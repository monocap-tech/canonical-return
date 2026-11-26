from pathlib import Path
import shutil,zipfile,hashlib,json
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';cb=w/'cr_int_t1d_cb';ca=w/'cr_int_t1d_ca_odd';out=r/'custody/ct/generated';out.mkdir(exist_ok=True)
# CK's exact BR production authentication; retained as a separate source step.
src=(r/'custody/ck/prepare.py').read_text().split("for n in ['bw_prepare.py'")[0].replace("custody/ck/generated","custody/ct/generated")
exec(compile(src,str(r/'custody/ck/prepare.py'),'exec'),dict(__file__=str(r/'custody/ck/prepare.py')))
bo=w/'cr_int_t1d_bo';bo.mkdir(exist_ok=True)
with zipfile.ZipFile(r/'ce-output/archives/CT-INT-T1D-BO-reproducibility.zip') as z:
    for name in ['odd_joint.npz','odd_overlap.npz']:
        hits=[n for n in z.namelist() if Path(n).name==name]
        assert len(hits)==1,hits
        (bo/name).write_bytes(z.read(hits[0]))
        print('BO source',hits[0],hashlib.sha256((bo/name).read_bytes()).hexdigest(),flush=True)
for name in ['bn_energy.py','bq_corrector.py','bq_tail.py','odd_merged_corrector.npz']:
    shutil.copyfile(ca/name,cb/name)
# Reuse the original sampled inverse solve, then select unused modes by the
# archived BY residual-overlap heuristic, targeted solely at direction 16.
for n in ['odd_merged_corrector.npz','odd_combined_projection.npz','odd_merged_tail.npz','odd_BQ_energy.npz']:
    shutil.copyfile(r/'custody/cr/generated/merged'/n,cb/n)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('d=119;m=64')==1
s=s.replace('d=119;m=64','d=119;m=128')
assert s.count('fac=1.02;')==1
s=s.replace('fac=1.02;','fac=1.0;')
s=s[:s.index('# represent X as full')]
s+="""
ev,Q=np.linalg.eigh(M);q=Q[:,-16:-15];a=np.vstack([q,X@q,-V@q]);u=U@a
old=np.load('odd_merged_corrector.npz')['indices']
cand=np.arange(3,4096,2);cand=cand[~np.isin(cand,old)];score=np.empty(len(cand))
for i in range(0,len(cand),128):
    ns=cand[i:i+128];basis=np.sin(np.pi*xx[:,None]*ns[None,:]/4)/np.sqrt(2)
    score[i:i+len(ns)]=np.linalg.norm(basis.T@u*8/N,axis=1)*np.log(2+ns)**.5
order=np.argsort(score)[-16:][::-1];new=cand[order]
assert len(old)==267 and len(new)==16 and not np.intersect1d(new,old).size
print('CT targets',ev[-16:-15].tolist(), 'selected',new.tolist(), 'scores',score[order].tolist(),flush=True)
np.save('odd_new.npy',new);np.savez_compressed('ct_selection.npz',cand=cand,score=score,eig=ev[-16:-15],new=new,direction=q)
"""
(cb/'ct_select.py').write_text(s)
shutil.copyfile(cb/'ct_select.py',out/'ct_select.py')
