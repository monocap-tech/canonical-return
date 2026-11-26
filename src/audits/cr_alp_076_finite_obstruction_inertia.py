#!/usr/bin/env python3
"""CR-ALP-076 finite scalar obstruction inertia.
Not an infinite Schur certificate. Finite residual Gram is a lower bound;
its scalar comparison is an UPPER bound on the complete scalar comparison.
No original Weil negativity follows.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES
S=10**18
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def interval(A,i,j):
    a,b=map(int,A[i][j]);c,d=map(int,A[j][i])
    l,h=max(a,c),min(b,d);assert l<=h
    return l,h
def bracket(w,use_in):
    K,G,I=w['K'],w['Gout'],w['Gin'];D=int(w['scale']);n=len(K)
    M=[];R=[]
    for i in range(n):
        row=[];err=0
        for j in range(n):
            kl,kh=interval(K,i,j);gl,gh=interval(G,i,j)
            if use_in:
                il,ih=interval(I,i,j);gl+=il;gh+=ih
            l=floor(F((51*kl-200*gh)*S,51*D))
            h=ceil(F((51*kh-200*gl)*S,51*D))
            assert l<=h
            m=(l+h)//2;row.append(m);err+=max(m-l,h-m)
        M.append(row);R.append(err)
    assert all(M[i][j]==M[j][i] for i in range(n) for j in range(n))
    lo=[r[:] for r in M];hi=[r[:] for r in M]
    for i in range(n):lo[i][i]-=R[i];hi[i][i]+=R[i]
    return lo,hi
def inertia(A):
    B=[r[:] for r in A];prev=1;pos=neg=0;negative_pivot_positions=[]
    for k in range(len(B)):
        p=B[k][k];assert p!=0,'Zero leading pivot; pivoting required'
        if p*prev>0:pos+=1
        else:neg+=1;negative_pivot_positions.append(k)
        for i in range(k+1,len(B)):
            for j in range(i,len(B)):
                q,r=divmod(p*B[i][j]-B[i][k]*B[k][j],prev)
                assert r==0;B[i][j]=B[j][i]=q
        prev=p
    return {'positive':pos,'negative':neg,'zero':0,
            'negative_elimination_pivot_positions':negative_pivot_positions}
def selftest():
    assert inertia([[2,1],[1,-2]])['negative']==1
    assert inertia([[-2,1],[1,-2]])['negative']==2
    assert inertia([[2,1],[1,2]])['positive']==2
def main(p,out):
    selftest();result={}
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        n=len(w['K']);assert n==(43 if k=='1' else 42)
        result[k]={}
        for use_in in (False,True):
            lo,hi=bracket(w,use_in)
            a,b=inertia(lo),inertia(hi)
            assert a['positive']<=b['positive'] and a['negative']>=b['negative']
            result[k]['Gin_plus_Gout' if use_in else 'Gout_only']={
              'lower':a,'upper':b,
              'exact_finite_inertia_determined':all(a[t]==b[t] for t in ('positive','negative','zero'))}
    r={'status':'FINITE_SCALAR_COMPARISON_ONLY','sectors':result,
       'scope':'Frozen 056 trial, K1250 finite residual. Negative elimination pivot positions are not physical coordinate witnesses. No new index cap or infinite Schur positivity.'}
    out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
