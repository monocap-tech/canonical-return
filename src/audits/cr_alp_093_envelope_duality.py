#!/usr/bin/env python3
"""CR-ALP-093 staged: two-constraint robust envelope duality.

For any PSD B with B<=M and Z^T B Z<=A, rank-one vv^T<=B
can reproduce x^T B x for any fixed x (v=Bx/sqrt(x^T B x)).
Thus robust B<K reduces to the ellipsoid intersection
  v^T M^-1 v<=1, (Z^T v)^T A^-1(Z^T v)<=1.
The real SDP relaxation of max v^T K^-1 v has two scalar
constraints; an extreme optimal solution has rank <=1.
A positive convex inverse metric certifies the robust bound.
A feasible rational rank-one counterexample disproves only this
relaxed envelope, not the actual inverse operator or Weil form.
Numerical optimization is candidate generation only.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import numpy as np
from scipy.optimize import minimize_scalar
from cr_alp_086_terminal_moment import data,full,exact_spd
from cr_alp_085_obstruction_block import mm,transpose,inv
from cr_alp_075_shifted_energy import S,floor
from cr_alp_082_inertia_receiver import inertia
def fmat(A):return np.array([[float(v) for v in r] for r in A])
def outer(v):return [[a*b for b in v] for a in v]
def sub(A,B):return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def cert(K,M,Z,A,w):
    assert 0<=w<1
    T=[[ (1-w)*a+w*b for a,b in zip(r,s)] for r,s in zip(inv(M),mm(Z,mm(inv(A),transpose(Z))))]
    R=sub(mm(K,mm(T,K)),K)
    n=len(K);C=[[floor(v*S) for v in r] for r in R]
    for i in range(n):C[i][i]-=n-1
    ii=inertia(C)
    return {'weight':str(w),'inertia':ii,'full_pass':ii==[n,0,0]}
def search(K,M,Z,A):
    n=len(K);k,m,z,a=map(fmat,(K,M,Z,A))
    L=np.linalg.cholesky(k);Mi=np.linalg.inv(m);P=z@np.linalg.solve(a,z.T)
    H=L.T@Mi@L;J=L.T@P@L
    def eig(w):
        return np.linalg.eigh((1-w)*H+w*J)
    r=minimize_scalar(lambda w:-eig(w)[0][0],bounds=(0,1-1e-9),method='bounded',options={'xatol':1e-12})
    candidates=[0.,float(r.x),max(0.,float(r.x)-1e-6),min(1-1e-9,float(r.x)+1e-6)]
    w=max(candidates,key=lambda w:eig(w)[0][0]);lam,V=eig(w)
    q=F(min(10**12-1,max(0,round(w*10**12))),10**12)
    result={'numeric_dual_margin':float(lam[0]),'dual_weight':str(q)}
    if lam[0]>1:
        result['certificate']=cert(K,M,Z,A,q)
        if result['certificate']['full_pass']:return result
    # Find a strict feasible rational rank-one witness, if available.
    # In an eigenspace of dimension >1, mix candidate eigenvectors.
    props=[V[:,j] for j in range(min(n,4))]
    for j in range(1,min(n,4)):
        for t in np.linspace(0,2*np.pi,65):props.append(np.cos(t)*V[:,0]+np.sin(t)*V[:,j])
    best=None
    for u in props:
        v=L@u
        scale=max(float(v@Mi@v),float(v@P@v))
        if scale<=0:continue
        v=v/np.sqrt(scale)*.999
        gain=float(v@np.linalg.solve(k,v))
        if best is None or gain>best[0]:best=(gain,v)
    if best is not None:
        gain,v=best;result['numeric_rank_one_gain']=gain
        if gain>1:
            vv=[F(round(float(t)*10**12),10**12) for t in v]
            c=[r[0] for r in mm(transpose(Z),[[t] for t in vv])]
            if exact_spd(sub(M,outer(vv))) and exact_spd(sub(A,outer(c))):
                ii=inertia(sub(K,outer(vv)))
                if ii[1]>0:result['relaxation_counterexample']={'v':list(map(str,vv)),'inertia_K_minus_vvT':ii}
    return result
def main(root,prior,block,previous,received,out):
    g=json.loads(previous.read_text());assert g['rank']==83
    r=json.loads(received.read_text());assert r['status']=='CONDITIONAL_CAPPED_BLOCK_RECEIVER'
    report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        n=len(J)
        assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=n-1
        B=[list(map(F,row)) for row in r['sectors'][k]['block_majorant']];assert exact_spd(B)
        GL=[[F(int(v),S) for v in row] for row in g['sectors'][k]['correction_gram_lower_grid']]
        GU=[[F(int(v),S) for v in row] for row in g['sectors'][k]['correction_gram_upper_grid']]
        K=[[F(H[i][j],S)+GL[i][j] for j in range(n)] for i in range(n)];assert exact_spd(K)
        BP=mm(transpose(W),mm(A,W))
        results=[];best=max(n-1,r['sectors'][k]['rank'])
        for t in [F(0),F(1,20),F(1,10),F(1,5),F(3,10)]:
            M=[[F(D[i],S)*(i==j)+GU[i][j] for j in range(n)] for i in range(n)] if not t else [[(1+t)*(F(D[i],S)*(i==j)-DP[i][j])+(1+1/t)*BP[i][j]+GU[i][j] for j in range(n)] for i in range(n)]
            assert exact_spd(M)
            res=search(K,M,Z,B);res['global_eta']=str(t);results.append(res)
            if res.get('certificate',{}).get('full_pass'):best=n;break
        report[k]={'rank':best,'tests':results}
        print(k,best,[(x['numeric_dual_margin'],bool(x.get('relaxation_counterexample'))) for x in results],flush=True)
    rank=sum(v['rank'] for v in report.values())
    out.write_text(json.dumps({'status':'CONDITIONAL_ENVELOPE_DUALITY','rank':rank,'conditional_index_upper':85-rank,'sectors':report,
    'premises':'Authenticated 091/092 bounds and inherited original source and continuum premises. Rejection of a two-constraint relaxation does not reject simultaneous additional constraints, actual operator positivity or RH.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
