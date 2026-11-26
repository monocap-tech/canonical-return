from pathlib import Path
import shutil,zipfile,hashlib,json
r=Path(__file__).resolve().parents[2];w=r/'ce-output/pristine';cb=w/'cr_int_t1d_cb';ca=w/'cr_int_t1d_ca_odd';out=r/'custody/dv/generated';out.mkdir(exist_ok=True)
# CK's exact BR production authentication; retained as a separate source step.
src=(r/'custody/ck/prepare.py').read_text().split("for n in ['bw_prepare.py'")[0].replace("custody/ck/generated","custody/dv/generated")
exec(compile(src,str(r/'custody/ck/prepare.py'),'exec'),dict(__file__=str(r/'custody/ck/prepare.py')))
bo=w/'cr_int_t1d_bo';bo.mkdir(exist_ok=True)
with zipfile.ZipFile(r/'ce-output/archives/CR-INT-T1D-BO-reproducibility.zip') as z:
    for name in ['odd_joint.npz','odd_overlap.npz']:
        hits=[n for n in z.namelist() if Path(n).name==name]
        assert len(hits)==1,hits
        (bo/name).write_bytes(z.read(hits[0]))
        print('BO source',hits[0],hashlib.sha256((bo/name).read_bytes()).hexdigest(),flush=True)
for name in ['bn_energy.py','bq_corrector.py','bq_tail.py','odd_merged_corrector.npz']:
    shutil.copyfile(ca/name,cb/name)
# Reuse the original sampled inverse solve, then select unused modes by the
# archived BY residual-overlap heuristic, jointly targeting directions 31 through 34.
for n in ['odd_merged_corrector.npz','odd_combined_projection.npz','odd_merged_tail.npz','odd_BQ_energy.npz']:
    shutil.copyfile(r/'custody/du/generated/merged'/n,cb/n)
p=cb/'bw_propose_full.py';s=p.read_text()
assert s.count('d=119;m=64')==1
s=s.replace('d=119;m=64','d=119;m=336')
assert s.count('fac=1.02;')==1
s=s.replace('fac=1.02;','fac=1.0;')
s=s[:s.index('# represent X as full')]
s+="""
ev,Q=np.linalg.eigh(M);print('DV prior spectrum top 36',ev[-36:].tolist(),flush=True);q=Q[:,-34:-30][:,::-1];a=np.vstack([q,X@q,-V@q]);u=U@a
old=np.load('odd_merged_corrector.npz')['indices']
cand=np.arange(3,4096,2);cand=cand[~np.isin(cand,old)];score=np.empty((len(cand),4))
for i in range(0,len(cand),128):
    ns=cand[i:i+128];basis=np.sin(np.pi*xx[:,None]*ns[None,:]/4)/np.sqrt(2)
    score[i:i+len(ns)]=np.abs(basis.T@u*8/N)*np.log(2+ns)[:,None]**.5
# Round-robin selection gives each residual direction four distinct modes.
assert np.all(np.isfinite(score)) and np.all(np.max(score,axis=0)>0)
rank=np.argsort(-score,axis=0,kind='stable')
order=[];owners=[]
for round_id in range(4):
    for col in range(4):
        k=next(int(k) for k in rank[:,col] if int(k) not in order)
        order.append(k);owners.append(31+col)
order=np.asarray(order);new=cand[order]
assert owners==[31,32,33,34]*4
assert len(old)==475 and len(new)==16 and not np.intersect1d(new,old).size
print('DV targets',ev[-34:-30][::-1].tolist(), 'selected',new.tolist(), 'owners',owners,'scores',score[order].tolist(),flush=True)
np.save('odd_new.npy',new);np.savez_compressed('dv_selection.npz',cand=cand,score=score,eig=ev[-34:-30][::-1],new=new,direction=q,owners=np.asarray(owners))
"""
(cb/'dv_select.py').write_text(s)
shutil.copyfile(cb/'dv_select.py',out/'dv_select.py')
