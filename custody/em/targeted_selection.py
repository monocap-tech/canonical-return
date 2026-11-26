"""EM targeted selection; numerical scores are not certificates."""
from pathlib import Path
import hashlib,json,shutil,zipfile
import numpy as np
root=Path(__file__).resolve().parents[2];w=root/'ce-output/pristine';cb=w/'cr_int_t1d_cb';ek=root/'custody/ek/generated';out=root/'custody/em/generated';out.mkdir(parents=True,exist_ok=True)
for line in (ek/'sha256.txt').read_text().splitlines():
    digest,name=line.split(maxsplit=1)
    with (ek/name.removeprefix('./')).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==digest
report=json.loads((ek/'screen.json').read_text())
assert report['ideal']['proxy32']<0 and report['best']['proxy32']<0
p=np.load(ek/'candidate32.npz');X,V,Q=[p[k] for k in ('X','V','Q')]
assert X.shape==(352,119) and Q.shape==(119,32) and np.max(abs(Q.T@Q-np.eye(32)))<1e-8
screen=root/'custody/ek/screen.py';s=screen.read_text();marker='L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float);U=np.empty_like(Z)'
assert s.count(marker)==1
env=dict(__file__=str(screen),__name__='em_selection_operator')
exec(compile(s.split(marker)[0],str(screen),'exec'),env)
N,xx,mask,Z,mv,LinearOperator,cg=[env[k] for k in ('N','xx','mask','Z','mv','LinearOperator','cg')]
assert Z.shape[1]==119+352+V.shape[0]
q=Q[:,-2:];rhs=Z@np.vstack([q,X@q,-V@q]);del Z
L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float);u=np.empty_like(rhs);res=[]
for k in range(2):
    v,info=cg(L,rhs[:,k],rtol=1e-10,maxiter=240);assert info==0
    res.append(float(np.linalg.norm(mv(v)-rhs[:,k])/max(np.linalg.norm(rhs[:,k]),1e-30)));u[:,k]=v
old=np.load(cb/'odd_merged_corrector.npz')['indices'];assert len(old)==491
cand=np.arange(3,16384,2);cand=cand[~np.isin(cand,old)]
full=np.zeros((N,2));full[mask]=u;ff=np.fft.fft(full,axis=0);x0=-4+4/N
overlap=-(ff[cand]*np.exp(-1j*np.pi*cand*x0/4)[:,None]).imag*8/(N*np.sqrt(2))
for k in (0,len(cand)//2,len(cand)-1):
    direct=np.sin(np.pi*xx*cand[k]/4)@u*8/(N*np.sqrt(2))
    assert np.max(abs(overlap[k]-direct))<1e-10
score=np.abs(overlap)*np.sqrt(np.log(2+cand))[:,None]
assert np.all(np.isfinite(score)) and np.all(np.max(score,axis=0)>0)
rank=np.argsort(-score,axis=0,kind='stable');chosen=[];owners=[]
for k in [1,1,1,0]*4:
    i=next(int(i) for i in rank[:,k] if int(i) not in chosen);chosen.append(i);owners.append(31+k)
new=cand[chosen];assert len(new)==len(set(new))==16 and not np.intersect1d(new,old).size
np.save(cb/'odd_new.npy',new);np.save(out/'odd_new.npy',new)
np.savez_compressed(out/'selection.npz',cand=cand,score=score,overlap=overlap,new=new,owners=np.array(owners),Q=q)
result=dict(scope='sampled overlap heuristic only',candidate_count=len(cand),previous_support=len(old),selected=new.tolist(),owners=owners,max_sampled_relative_residual=max(res),selected_scores=score[chosen].tolist(),max_score=score.max(axis=0).tolist(),high_band_selected=int(np.count_nonzero(new>=4096)),old_proxy32=report['best']['proxy32'])
(out/'selection.json').write_text(json.dumps(result,indent=2)+'\n')
ca=w/'cr_int_t1d_ca_odd';bo=w/'cr_int_t1d_bo'
with zipfile.ZipFile(root/'ce-output/archives/CR-INT-T1D-BO-reproducibility.zip') as z:
    for name in ('odd_joint.npz','odd_overlap.npz'):
        hits=[n for n in z.namelist() if Path(n).name==name];assert len(hits)==1
        (bo/name).write_bytes(z.read(hits[0]))
for name in ('bn_energy.py','bq_corrector.py','bq_tail.py'):shutil.copyfile(ca/name,cb/name)
print('EM SELECTION',result,flush=True)
