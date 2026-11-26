#!/usr/bin/env python3
"""CR-ALP-075: exact shifted energy floors and finite residual rejection tests.

Source-only computation. K is the pinned corrected energy matrix;
Gout is a LOWER bound on the complete residual Gram. Never use Gout
as a global upper bound. No new positivity or index claim.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_075_nested_selection import select
S=10**18
GRID=10000
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def positive(A):
    B=[r[:] for r in A];prev=1
    for k in range(len(B)):
        p=B[k][k]
        if p<=0:return False
        for i in range(k+1,len(B)):
            for j in range(i,len(B)):
                num=p*B[i][j]-B[i][k]*B[k][j]
                q,r=divmod(num,prev);assert r==0
                B[i][j]=B[j][i]=q
        prev=p
    return True
def lower_matrix(w,J):
    K=w['K'];D=int(w['scale']);n=len(J)
    M=[];R=[]
    for i in J:
        row=[];radius=[]
        for j in J:
            a,b=map(int,K[i][j]);c,d=map(int,K[j][i])
            a,b=max(a,c),min(b,d);assert a<=b
            lo,hi=floor(F(a*S,D)),ceil(F(b*S,D))
            m=floor(F((a+b)*S,2*D))
            row.append(m);radius.append(max(m-lo,hi-m))
        M.append(row);R.append(sum(radius))
    assert all(M[i][j]==M[j][i] for i in range(n) for j in range(n))
    for i in range(n):M[i][i]-=R[i]
    return M
def test(M,q):
    assert S%GRID==0
    return positive([[v-(q*S//GRID if i==j else 0)
                      for j,v in enumerate(row)] for i,row in enumerate(M)])
def shift(M):
    assert test(M,0)
    l,h=0,GRID+1
    while h-l>1:
        mid=(l+h)//2
        if test(M,mid):l=mid
        else:h=mid
    assert test(M,l)
    return F(l,GRID)
def run(p,out):
    result={}
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        sel=select(w,k);J=sel['complete_order'];D=int(w['scale'])
        G=w['Gout'];rows={}
        for n in range(10,len(J)+1):
            I=J[:n];M=lower_matrix(w,I)
            lam=shift(M)
            trace=sum(max(F(0),F(int(G[i][i][0]),D)) for i in I)
            budget=ELL*lam
            rows[str(n)]={'shifted_K_lower':str(lam),
                          'finite_residual_trace_lower':str(trace),
                          'global_residual_sum_budget':str(budget),
                          'finite_trace_already_exceeds_this_budget':trace>=budget}
        result[k]={'coordinate_order':J,'prefixes':rows}
    r={'status':'SOURCE_ENERGY_ONLY','coercivity':str(ELL),
       'sectors':result,
       'premises':'Pinned corrected energy and finite projected residual matrices. Finite trace is a lower bound only. No new physical integrations, index reduction or Lean proof.'}
    out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))
