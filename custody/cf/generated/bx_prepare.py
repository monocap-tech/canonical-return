from bu_circle import *
from math import isqrt
bk=parent/'cr_int_t1d_bk'/'even';bp=parent/'cr_int_t1d_bp';bu=parent/'cr_int_t1d_bu';a=np.load(root/'bx_full_proposal.npz');Q=a['Q'][:,:4];X=a['X'];t=a['t'][:,:4];v0=a['v0'][:,:4];g0=a['g0'][:,:4];B=np.load(root/'even_combined_projection.npz');Y=B['Y'];x=X.astype(object)@Q.astype(object);Bq=np.vstack([Q.astype(object)*S,x]);yq=Y.astype(object)@Bq;yn=(yq//(S*S)).astype(np.int64)
physical=np.load(bk/'physical_dictionary.npz');E=np.load(bk/'exceptional.npz');ES=int(E['scale']);j=physical['indices'];rj=a['rj'];alpha=a['alpha'];assert not np.intersect1d(j,rj).size
wroot=np.array([isqrt(int(v)*S*S//int(physical['weight_scale'])) for v in physical['weights']],object);aroot=np.array([isqrt(int(v)*S*S//10**8) for v in alpha],object);raw=E['E'].astype(object)@v0.astype(object);fc=-(wroot[:,None]*raw//(ES*S));rc=-(aroot[:,None]*t.astype(object)//S);gc=-g0
ms=np.load(bp/'even_head_full.npz')['indices'];cy=yn.astype(object)*isqrt(2*S*S)*((-1)**(ms//2))[:,None]//S;ids=np.union1d(np.r_[j,rj],2*ms);C=np.zeros((len(ids),Q.shape[1]),object);C[np.searchsorted(ids,j)]+=fc;C[np.searchsorted(ids,rj)]+=rc;C[np.searchsorted(ids,2*ms)]+=cy;assert max(abs(x) for x in C.flat)<2**63
errs=[]
for k in range(Q.shape[1]):
 eyr=norm(yq[:,k:k+1]-yn[:,k:k+1].astype(object)*S*S,S**3);ef=norm(raw[:,k:k+1],ES*S)/S+F(len(j),S);er=norm(t[:,k:k+1],S)/S+F(len(rj),S);ec=norm(yn[:,k:k+1],S)/S+F(len(ms),S);errs.append(eyr+ef+er+ec)
print('representation PASS',len(ids),'max error',float(max(errs)),flush=True)
B0=np.load(bp/'even_projected.npz');np.savez_compressed(root/'bu_trial.npz',Q=Q,X=X,t=t,alpha=alpha,js=ids,c=C.astype(np.int64),gc=gc,errors=np.array([str(e) for e in errs]),scale=S,eps=str(B0['eps']),mu=str(B0['mu']),tail=str(B0['tail']),head_error=str(B['head_error']),Uerror=str(B['U_error']),Utail=str(np.load(root/'even_merged_tail.npz')['tail'][list(np.load(root/'even_merged_tail.npz')['cutoff']).index(16384)]));np.savez_compressed(root/'even_joint_proposal.npz',j=a['j'],T=a['T'][:,:4],scale=S)
