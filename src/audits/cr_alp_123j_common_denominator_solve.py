#!/usr/bin/env python3
"""Exact common-determinant solve for candidate recovery.

No statement about positive definiteness of the original interval matrix.
All candidates must still pass the original interval residual/energy receiver.
"""
from fractions import Fraction as F

def solve_common(A,b):
    n=len(A);assert n and len(b)==n and all(len(r)==n for r in A)
    assert A==list(map(list,zip(*A)))
    B=[list(r)+[v] for r,v in zip(A,b)];p0=1
    for k in range(n-1):
        p=B[k][k];assert p>0
        for i in range(k+1,n):
            for j in range(k+1,n+1):
                v,r=divmod(p*B[i][j]-B[i][k]*B[k][j],p0)
                assert r==0;B[i][j]=v
            B[i][k]=0
        p0=p
    den=B[-1][-2];assert den>0
    x=[0]*n
    for i in reversed(range(n)):
        x[i],r=divmod(den*B[i][-1]-sum(B[i][j]*x[j] for j in range(i+1,n)),B[i][i])
        assert r==0
    assert all(sum(a*v for a,v in zip(row,x))==den*t for row,t in zip(A,b))
    return x,den

def round_quotient(n,d):
    assert d>0
    q,r=divmod(n,d)
    return q+int(2*r>d or (2*r==d and q%2==1))

def trial(A,b,D=10**90):
    x,den=solve_common(A,b)
    return [F(round_quotient(v*D,den),D) for v in x]

def reduced_midpoint_trial(A,b,source_scale,target_scale=10**110,D=10**90):
    # Proposal only: a coarser exact rational system. Its minimizer need
    # not equal the original midpoint minimizer. Never use its energy in
    # place of the original outward interval evaluation.
    assert source_scale>=target_scale>0
    P=[[round_quotient(v*target_scale,source_scale) for v in r] for r in A]
    rhs=[round_quotient(v*target_scale,source_scale) for v in b]
    return trial(P,rhs,D)

if __name__=='__main__':
    from random import Random
    from cr_alp_123b_pole_residual import solve
    rng=Random(20261126)
    for n,scale in [(3,1),(12,10**30),(32,10**55),(48,10**55)]:
        X=[[rng.randrange(-scale,scale+1) for _ in range(n)] for _ in range(n)]
        A=[[sum(X[k][i]*X[k][j] for k in range(n))+(scale if i==j else 0) for j in range(n)] for i in range(n)]
        b=[rng.randrange(-scale,scale+1) for _ in range(n)]
        expected=[F(round(v*10**90),10**90) for v in solve(A,b)]
        assert trial(A,b)==expected
        print('PASS exact rounded coordinates',n,flush=True)
