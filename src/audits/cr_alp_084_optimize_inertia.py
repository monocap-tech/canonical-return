#!/usr/bin/env python3
"""CR-ALP-084: optimize full inertia, not coordinate-prefix positivity.
Floating proposals only. Every reported dimension is checked by exact
rational diagonal budget and symmetric congruence.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import numpy as np
from cr_alp_074_balanced_selection import checked,HASHES
from cr_alp_075_shifted_energy import lower_matrix,S,ceil
from cr_alp_083_frozen_inertia import inertia

def search(H,aa,initial,rank):
    n=len(H);A=np.asarray(H,dtype=float)/S
    a=np.asarray([float(x) for x in aa])
    assert np.all(a>0)
    rng=np.random.default_rng(8401+n)
    def weights(d):
        w=a/np.array([float(x)/S for x in d])
        return w/w.sum()
    def exact(w):
        # Integer weights normalized exactly, no floating-point acceptance.
        W=[max(1,int(x*10**14)) for x in w];T=sum(W)
        D=[max(1,ceil(x*S*T/F(z))) for x,z in zip(aa,W)]
        theta=sum(x*S/F(d) for x,d in zip(aa,D))
        assert theta<=1
        ii=inertia([[x-(D[i] if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(H)])
        return {'inertia':ii,'diagonal_grid':list(map(str,D)),'theta_upper':str(theta)}
    def eig(w,k):
        values,V=np.linalg.eigh(A-np.diag(a/w))
        return values[k],V[:,k]
    w0=weights(initial)
    base={'inertia':inertia([[x-(initial[i] if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(H)]),
          'diagonal_grid':list(map(str,initial)),
          'theta_upper':str(sum(x*S/F(d) for x,d in zip(aa,initial)))}
    assert F(base['theta_upper'])<=1 and base['inertia'][0]==rank
    best=base;diagnostics=[]
    # For target rank r, maximize eigenvalue n-r of the conservative
    # comparison. Softmax parameterization enforces positive weights.
    def soft(z):
        z=z-np.max(z);w=np.exp(z);return w/w.sum()
    for target in range(rank+1,n+1):
        k=n-target;proposals=[]
        starts=[weights(list(map(int,best['diagonal_grid']))),w0,np.sqrt(a)/np.sqrt(a).sum(),a/a.sum()]
        for j in range(12):
            z=np.log(w0)+rng.normal(0,.25 if j<6 else .7,n)
            starts.append(soft(z))
        for start in starts:
            z=np.log(start);val=-float('inf')
            for it in range(350):
                w=soft(z);v,x=eig(w,k)
                if v>val:val=v;bestw=w.copy()
                # d lambda / d w_i = a_i*x_i^2/w_i^2; chain through softmax.
                g=a*x*x/w;g=g-w*g.sum()
                norm=np.linalg.norm(g)
                if norm<1e-10:break
                direction=g/norm;t=.5
                while t>1e-7:
                    zz=np.clip(z+t*direction,-40,40)
                    nv,_=eig(soft(zz),k)
                    if nv>v+1e-5*t*norm:z=zz;break
                    t*=.5
                if t<=1e-7:break
            proposals.append((val,bestw))
        proposals.sort(key=lambda x:x[0],reverse=True)
        diagnostic={'target':target,'best_numeric_eigenvalue':float(proposals[0][0]),'exact_trials':[]}
        for val,w in proposals[:4]:
            row=exact(w);diagnostic['exact_trials'].append(row)
            if row['inertia'][0]>best['inertia'][0]:best=row
        diagnostics.append(diagnostic)
        if best['inertia'][0]<target:break
    return best,diagnostics

def main(root,prior,frozen,out):
    P=json.loads(prior.read_text());Q=json.loads(frozen.read_text())
    assert P['rank']==72 and Q['rank']==79
    res={}
    for k,h in HASHES.items():
        x=P['sectors'][k];y=Q['sectors'][k];J=x['coordinate_order']
        assert J==y['coordinate_order']
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',h)
        H=lower_matrix(w,J);aa=[F(x['inverse_bounds'][str(j)]) for j in J]
        initial=list(map(int,y['diagonal_grid']))
        best,diag=search(H,aa,initial,y['inertia'][0])
        res[k]={'best':best,'diagnostics':diag,'coordinate_order':J,'scale':str(S)}
    rank=sum(x['best']['inertia'][0] for x in res.values())
    result={'status':'CONDITIONAL_OPTIMIZED_INERTIA','rank':rank,'conditional_index_upper':85-rank,'sectors':res,
         'premises':'Same authenticated original fixed-collar matrix and CR-ALP-081 inverse bounds as 083. Floating search only; exact rational acceptance. No unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
