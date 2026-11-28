#!/usr/bin/env python3
"""CR-ALP-074 exact full-energy fallback; no extra physical integration.

Given the authenticated K bank and norm receiver output, test
max_{s in {+-1}^n} (s*b)^T H^-1 (s*b) < ell,
where H <= K is an exact symmetric rational lower matrix.
For every x, ||R x|| <= sum b_i |x_i| = max_s (s*b)^T x.
Cauchy-Schwarz in H gives R*L^-1 R <= (max/ell) H.
The logged b_i must already enclose actual projected residual norms.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES
ELL=F(51,200)
S=10**18
def floor(x):return F(x*S//1,S)
def ceil(x):return -floor(-x)
def lower_matrix(w,J):
    K=w['K'];scale=int(w['scale']);n=len(J)
    H=[[F(0) for j in J] for i in J];r=[F(0) for i in J]
    for i,a in enumerate(J):
        for j,b in enumerate(J):
            l,h=map(int,K[a][b]);ll,hh=map(int,K[b][a])
            l,h=max(l,ll),min(h,hh);assert l<=h
            lo,hi=floor(F(l,scale)),ceil(F(h,scale))
            H[i][j]=floor((lo+hi)/2)
            r[i]+=max(H[i][j]-lo,hi-H[i][j])
    assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
    for i in range(n):H[i][i]-=r[i]
    return H
def ldl(H):
    n=len(H);L=[[F(int(i==j)) for j in range(n)] for i in range(n)];D=[]
    for i in range(n):
        for j in range(i):
            L[i][j]=(H[i][j]-sum(L[i][k]*D[k]*L[j][k] for k in range(j)))/D[j]
        d=H[i][i]-sum(L[i][k]**2*D[k] for k in range(i))
        assert d>0,'Energy lower matrix not SPD'
        D.append(d)
    assert all(H[i][j]==sum(L[i][k]*D[k]*L[j][k] for k in range(min(i,j)+1)) for i in range(n) for j in range(n))
    return L,D
def vertex_max(H,b):
    n=len(b);L,D=ldl(H)
    assert all(t>=0 for t in b)
    # Transform each column of diag(b) through the unit lower factor.
    # Round the resulting quadratic matrix outward before enumeration,
    # retaining exact comparison with a 10^-18 grid.
    Y=[[F(0) for j in range(n)] for i in range(n)]
    for j in range(n):
        for i in range(n):
            Y[i][j]=(b[j] if i==j else F(0))-sum(L[i][k]*Y[k][j] for k in range(i))
    T=[[sum(Y[k][i]*Y[k][j]/D[k] for k in range(n)) for j in range(n)] for i in range(n)]
    diag=sum(ceil(T[i][i]) for i in range(n))
    pairs=[(i,j,floor(2*T[i][j]),ceil(2*T[i][j])) for i in range(n) for j in range(i)]
    best=None;sign=None
    for tail in product((-1,1),repeat=n-1):
        s=(1,)+tail
        upper=diag+sum(hi if s[i]==s[j] else -lo for i,j,lo,hi in pairs)
        if best is None or upper>best:best,sign=upper,s
    assert best is not None and best>=0
    return best,sign,min(D)
def selftest():
    H=[[F(2),F(1)],[F(1),F(2)]]
    m,s,d=vertex_max(H,[F(1,10),F(1,5)])
    assert F(7,150)<=m<=F(7,150)+F(4,S)
    assert s==(1,-1)
    m,_,_=vertex_max([[F(2),F(0)],[F(0),F(3)]],[F(1,10),F(1,5)])
    assert F(11,600)<=m<=F(11,600)+F(4,S)
def run(p,report,out):
    selftest();data=json.loads(report.read_text());result={}
    assert data['status'] in ('PASS_CONDITIONAL_RANK20','UNDECIDED')
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        rows=data['sectors'][k]['columns']
        J=[r['coordinate_zero_based'] for r in rows];assert len(set(J))==len(J)==10
        b=[F(r['actual_residual_norm_upper']) for r in rows]
        H=lower_matrix(w,J)
        m,s,d=vertex_max(H,b)
        trace=F(data['sectors'][k]['residual_norm_squared_sum_upper'])/F(data['sectors'][k]['K_lower'])
        assert sum(q*q for q in b)==F(data['sectors'][k]['residual_norm_squared_sum_upper'])
        # Either independent sufficient comparison can be used.
        budget=min(m,trace);gap=ELL-budget
        result[k]={'rank':10,'coordinates':J,'energy_pivot_lower':str(d),
            'sign_vertex_upper':str(m),'maximizing_sign_upper':s,
            'previous_trace_upper':str(trace),'selected_upper':str(budget),
            'coercivity_margin':str(gap),'status':'PASS_CONDITIONAL' if gap>0 else 'UNDECIDED'}
    ok=all(r['status']=='PASS_CONDITIONAL' for r in result.values())
    r={'status':'PASS_CONDITIONAL_RANK20' if ok else 'UNDECIDED',
       'sectors':result,'conditional_index_cap_if_received':65 if ok else 66,
       'premises':'Same authenticated source and actual norm bounds as 074 norm receiver. No extra integration or independent continuum proof.'}
    out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':
    if sys.argv[1:]==['--selftest']:selftest();print('PASS exact synthetic vertex tests')
    else:run(*(Path(x) for x in sys.argv[1:]))
