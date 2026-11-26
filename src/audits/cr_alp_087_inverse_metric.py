#!/usr/bin/env python3
"""CR-ALP-087: combine PSD inverse-energy constraints in the inverse metric.

For B>=0, U* B U <= A with A>0 implies
B^(1/2) U A^-1 U* B^(1/2) <= I (equal nonzero singular values).
Any convex combination T of U A^-1 U* has B^(1/2) T B^(1/2)<=I.
If T>0, B<=T^-1, even for singular B. Thus K_lower>T^-1
suffices. Inertia(K_lower-T^-1)=inertia(K_lower T K_lower-K_lower)
for K_lower,T>0 by the two Schur complements of [[K_lower,I],[I,T]].

Cutting-plane LP is only a floating witness search. Exact rational
weights, reconstructed moment bounds and integer congruence receive it.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import numpy as np
from scipy.optimize import linprog
from cr_alp_086_terminal_moment import data,full,entrymap,input_for,exact_spd
from cr_alp_085_obstruction_block import mm,transpose,inv
from cr_alp_080_selected_second_moment import receive as moments
from cr_alp_078_moment_receiver import sha
from cr_alp_075_shifted_energy import S,floor,positive
from cr_alp_082_inertia_receiver import inertia

def fmat(A):return np.array([[float(x) for x in r] for r in A])
def outer(a,b):return [[x*y for y in b] for x in a]
def search(H,terms):
    h=fmat(H)/S;n=len(H)
    L=np.linalg.cholesky(h)
    M=np.array([L.T@fmat(x)@L for x in terms]);m=len(M)
    # All weights nonnegative, first constraint has strictly positive
    # weight so T is SPD. Numerical LP only proposes convex combinations.
    bounds=[(1e-8,1)]+[(0,1)]*(m-1)+[(None,None)]
    equal=np.array([[1.]*m+[0.]])
    obj=np.array([0.]*m+[-1.])
    cuts=[]
    best=[]
    def evaluate(w):
        ev,V=np.linalg.eigh(np.einsum('i,ijk->jk',w,M))
        best.append((float(ev[0]),w.copy()))
        for j in range(min(4,n)):
            v=V[:,j];cuts.append(np.r_[-np.einsum('i,kil,l->k',v,M,v),1.])
        return float(ev[0])
    w=np.zeros(m);w[0]=1;evaluate(w)
    w=np.ones(m)/m;evaluate(w)
    last=None
    for it in range(300):
        r=linprog(obj,A_ub=np.array(cuts),b_ub=np.zeros(len(cuts)),
          A_eq=equal,b_eq=[1.],bounds=bounds,method='highs')
        if not r.success:break
        w=np.maximum(0,r.x[:-1]);w[0]=max(w[0],1e-8);w/=w.sum()
        val=evaluate(w);last=float(r.x[-1])
        if val>1.000001 or last-max(x[0] for x in best)<1e-7:break
    best.sort(key=lambda x:x[0],reverse=True)
    return best[:8],last,len(cuts)

def certify(H,terms,w):
    n=len(H);W=[max(0,int(round(x*10**12))) for x in w];W[0]=max(W[0],1)
    den=sum(W);weights=[F(x,den) for x in W]
    assert sum(weights)==1 and weights[0]>0 and min(weights)>=0
    T=[[sum(t*A[i][j] for t,A in zip(weights,terms) if t) for j in range(n)] for i in range(n)]
    # T>0 is automatic from the strictly positive full D^-1 weight.
    h=[[F(x,S) for x in row] for row in H]
    R=mm(h,mm(T,h))
    C=[[floor((R[i][j]-h[i][j])*S) for j in range(n)] for i in range(n)]
    assert all(C[i][j]==C[j][i] for i in range(n) for j in range(n))
    for i in range(n):C[i][i]-=n-1
    ii=inertia(C)
    return {'weights':list(map(str,weights)),'inertia':ii,'full_pass':ii==[n,0,0]}

def main(root,bank,prior,block,terminal,out):
    old=json.loads(prior.read_text());assert old['rank']==79
    p=json.loads((terminal/'result.json').read_text());assert p['rank']==83
    sel=json.loads((terminal/'terminal_selection.json').read_text())
    E=entrymap(bank);report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        n=len(J);assert positive(H)
        assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=n-1
        s=sel[k];q=list(map(F,s['q']))
        raw,eps=input_for(bank,E,k,J,Z,q);f=terminal/f'input_{k}_terminal.txt'
        assert raw==f.read_text() and sha(f)==s['input_sha256'] and eps==F(s['error_upper'])
        r=moments(terminal/f'mode_{k}_terminal.log',terminal/f'second_{k}_terminal.log',eps)
        assert r==p['sectors'][k]['moments']
        z=[sum(x*y for x,y in zip(row,q)) for row in Z]
        a=F(p['sectors'][k]['retained_inverse_upper'])
        base=sum(x*y for x,y in zip(q,[r[0] for r in mm(A,[[x] for x in q])]))
        assert a==min(base,F(r['inverse_upper'])) and a>0
        assert exact_spd(A)
        aa=[F(old['sectors'][k]['inverse_bounds'][str(j)]) for j in J]
        assert all(x>0 for x in aa) and sum(x*S/F(d) for x,d in zip(aa,D))<=1
        terms=[[[F(S,D[i]) if i==j else F(0) for j in range(n)] for i in range(n)],
               mm(Z,mm(inv(A),transpose(Z))),
               [[x/a for x in row] for row in outer(z,z)]]
        labels=['full_diagonal','physical_affine_block','terminal_quadratic']
        for k0,x in enumerate(aa):
            terms.append([[1/x if i==j==k0 else F(0) for j in range(n)] for i in range(n)])
            labels.append('coordinate_'+str(J[k0]))
        proposals,lp,cuts=search(H,terms)
        best=n-1;trials=[]
        for val,w in proposals:
            row=certify(H,terms,w);row['numeric_min_eigenvalue']=val;trials.append(row)
            best=max(best,row['inertia'][0])
            if best==n:break
        report[k]={'rank':best,'constraint_labels':labels,'trials':trials,
            'numeric_lp_upper_diagnostic':lp,'number_of_rayleigh_cuts':cuts,
            'terminal_inverse_upper':str(a),'pass_full_sector':best==n}
    rank=sum(x['rank'] for x in report.values())
    ans={'status':'CONDITIONAL_INVERSE_METRIC_RECEIVER','rank':rank,
      'conditional_index_upper':85-rank,'sectors':report,
      'premises':'Inherited continuum, source and numerical enclosure premises, authenticated CR-ALP-082/085/086 bounds; convex weights and congruence verified exactly. LP nonacceptance is not a mathematical obstruction or original Weil negative witness.'}
    out.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
