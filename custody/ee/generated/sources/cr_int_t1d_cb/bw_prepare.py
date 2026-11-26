from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np
root=Path(__file__).resolve().parent;parent=root.parent;bp=parent/'cr_int_t1d_bp';bk=parent/'cr_int_t1d_bk'/'odd';S=10**12
src=(bp/'bp_project.py').read_text();exec(src[src.index('obj=lambda'):src.index('h=np.load(root/')]);a=np.load(root/'bw_full_proposal.npz');Q=a['Q'];X=a['X'];t=a['t'];v0=a['v0'];g0=a['g0'];B=np.load(root/'odd_combined_projection.npz');Y=B['Y'];x=X.astype(object)@Q.astype(object);Bq=np.vstack([Q.astype(object)*S,x]);yq=Y.astype(object)@Bq;yn=(yq//(S*S)).astype(np.int64)
physical=np.load(bk/'physical_dictionary.npz');E=np.load(bk/'exceptional.npz');ES=int(E['scale']);j=physical['indices'];C=np.load(parent/'cr_int_t1d_br'/'odd_surplus_cross.npz');rj=C['indices'];assert np.array_equal(C['weights'],a['alpha']);assert not np.intersect1d(j,rj).size
wroot=np.array([isqrt(int(v)*S*S//int(physical['weight_scale'])) for v in physical['weights']],object);aroot=np.array([isqrt(int(v)*S*S//10**8) for v in C['weights']],object)
raw=E['E'].astype(object)@v0.astype(object);fc=-(wroot[:,None]*raw//(ES*S));rc=-(aroot[:,None]*t.astype(object)//S);c=np.vstack([fc,rc]).astype(np.int64);js=np.r_[j,rj];gc=-g0
# v0 denotes arbitrary projection on raw TN, not projected F. TN is in span(F,g), g arbitrary exact.
errs=[]
for k in range(Q.shape[1]):
 eyr=norm(yq[:,k:k+1]-yn[:,k:k+1].astype(object)*S*S,S**3)
 ef=norm(raw[:,k:k+1],ES*S)/S+F(len(j),S);er=norm(t[:,k:k+1],S)/S+F(len(rj),S);errs.append(eyr+ef+er)
print('representation max error',float(max(errs)),flush=True)
B0=np.load(bp/'odd_projected.npz');np.savez_compressed(root/'bs_trial.npz',Q=Q,X=X,t=t,yn=yn,v0=v0,g0=g0,js=js,c=c,gc=gc,errors=np.array([str(e) for e in errs]),scale=S,ms=np.load(bp/'odd_head_full.npz')['indices'],a=C['weights'],eps=str(B0['eps']),mu=str(B0['mu']))
np.savez_compressed(root/'odd_joint_proposal.npz',j=a['j'],T=a['T'],scale=S)
