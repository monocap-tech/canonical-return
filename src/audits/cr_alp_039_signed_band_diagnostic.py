#!/usr/bin/env python3
"""NUMERICAL ONLY: signed original 12-mode band coupling and Schur test.

The finite band's exact pole-annihilating projector gives a lower bound
on the *full* high-tail residual. These float results are NOT intervals
or proof of signs. The physical Weil Hmid input comes from CR-ALP-034.
"""
import sys,time
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'cr_alp_038'))
from combined_residual_pilot import H,M,cross_interior
from anchor_residual_pilot import anchor_coeff
from scipy.linalg import eigh
N=15000; MAX=45000; CHUNK=800; REM=3000

for block,ids in [('odd',[3,5,7,9,11,1]),('even',[4,6,8,10,12,2])]:
    j=np.arange(N+1+(block=='even'),MAX+1,2)
    columns=[]; cc=None
    start=time.time()
    for k in ids[:-1]:
        p,g,c=cross_interior(k,j)
        columns.append(p+g)
        cc=c
    ac=[];c_anchor=[]
    for sub in np.array_split(j,(len(j)+CHUNK-1)//CHUNK):
        all_,p,g,c,gamma=anchor_coeff(ids[-1],sub,rem_count=REM)
        ac.append(all_);c_anchor.append(c)
    aa=np.concatenate(ac)
    assert np.allclose(cc,np.concatenate(c_anchor),atol=1e-10,rtol=1e-10)
    columns.append(aa)
    T=np.stack(columns,axis=1)
    c=cc
    band = T.T@T-np.outer(T.T@c,T.T@c)/(c@c)
    band=(band+band.T)/2
    h=np.array([[H[i,k] for k in ids] for i in ids]);h=(h+h.T)/2
    m=np.array([[M[i,k] for k in ids] for i in ids]);m=(m+m.T)/2
    res_eig=np.linalg.eigvalsh(band)
    ev0=eigh(h,m,eigvals_only=True)[0]
    sch=h-4*band
    evs=eigh(sch,m,eigvals_only=True)
    vals,vecs=eigh(sch,m)
    u=vecs[:,0];g=u@h@u;b=u@band@u
    print(block,'jcount',len(j),'band anchor residual lower',np.sqrt(max(0,band[-1,-1])), 'band interior coupled norms',np.sqrt(np.maximum(0,np.diag(band)[:-1])))
    print('  H minimum generalized',ev0,'H-4 Bband generalized eigs',evs)
    print('  candidate negative min: H',g,'4Band',4*b,'gap',g-4*b)
    print('  vector:',dict(zip(ids,np.round(u,9))))
    print('  elapsed',time.time()-start)
    import json
    out={'block':block,'j_count':len(j),'modes':ids,'band_anchor_norm':float(np.sqrt(max(0,band[-1,-1]))),'eigen_H_minus_4_Bband':evs.tolist(),'gap_witness':float(g-4*b),'candidate':u.tolist(),'band_matrix':band.tolist()}
    (Path(__file__).resolve().parent/f'CR-ALP-039_{block}_band_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print('NUMERICAL PILOT ONLY. No interval certificate and no negative Weil source claimed.')
