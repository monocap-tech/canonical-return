#!/usr/bin/env python3
"""Conditioned block candidate generator (proposal only).
Input matrices are outward intervals of the original, not K or N.
Return exact rational coordinates; caller must certify original joint energy,
complete residual and conditioned response separately.
"""
from fractions import Fraction as F
import numpy as np
import mpmath as mp

def candidate(H,R,C,e,s,V,S,I,ZERO):
    n=len(H);m=len(C)
    assert len(R)==m==len(e)
    assert len(s)==len(V)==n
    assert all(len(v)==n for v in H+R+V)
    assert all(len(v)==m for v in C)
    # Transform first at original interval precision; never cast the
    # ill-conditioned raw columns or V to binary64 before cancellation.
    HV=[[sum((H[i][j]*V[j][k] for j in range(n)),ZERO)
         for k in range(n)] for i in range(n)]
    HH=[[sum((V[i][j]*HV[i][k] for i in range(n)),ZERO)
         for k in range(n)] for j in range(n)]
    RR=[[sum((r[i]*V[i][k] for i in range(n)),ZERO)
         for k in range(n)] for r in R]
    ss=[sum((V[i][j]*s[i] for i in range(n)),ZERO) for j in range(n)]
    def mid(v):return F(v.l+v.h,2*S)
    def m64(v):return float(mid(v))
    B=np.array([[m64(v) for v in row] for row in C])
    Q=np.array([[m64(v) for v in row] for row in RR])
    A=np.array([[m64(v) for v in row] for row in HH])
    E=np.array([m64(v) for v in e]);s0=np.array([m64(v) for v in ss])
    B=(B+B.T)/2;A=(A+A.T)/2
    L=np.linalg.cholesky(B)
    # Use one factorization, with multiple RHS.
    from scipy.linalg import cho_solve
    X=cho_solve((L,True),np.column_stack([Q,E]))
    T,x=X[:,:n],X[:,n]
    Schur=A-Q.T@T;Schur=(Schur+Schur.T)/2
    z=np.linalg.solve(Schur,-s0+Q.T@x)
    v=-x-T@z
    assert np.all(np.isfinite(z)) and np.all(np.isfinite(v))
    D=10**90
    zr=[F(round(F(float(x))*D),D) for x in z]
    vr=[F(round(F(float(x))*D),D) for x in v]
    # Freeze the midpoint transformation itself, so raw coordinate change
    # is an exact rational vector, not an interval-valued witness.
    raw=[sum((mid(V[i][j])*zr[j] for j in range(n)),F(0)) for i in range(n)]
    return raw,vr,zr

def selftest():
    # Small rational test exercises the actual interval-conditioned path
    # with a severely scaled raw coordinate and nontrivial mixed block.
    from cr_alp_046_full_finite import I,S,ZERO
    H=[[I(2),ZERO],[ZERO,I(F(3,10**40))]]
    V=[[I(1),ZERO],[ZERO,I(10**20)]]
    R=[[I(F(1,3)),I(F(1,10**21))]]
    C=[[I(4)]];e=[I(F(1,7))];s=[I(F(1,5)),I(F(1,10**20))]
    z,v,_=candidate(H,R,C,e,s,V,S,I,ZERO)
    low=[s[i]+sum((H[i][j]*I(z[j]) for j in range(2)),ZERO)+R[0][i]*I(v[0]) for i in range(2)]
    high=e[0]+sum((R[0][j]*I(z[j]) for j in range(2)),ZERO)+C[0][0]*I(v[0])
    cond=[sum((V[i][j]*low[i] for i in range(2)),ZERO) for j in range(2)]
    assert max(x.up() for x in cond+[high])<F(1,10**13)
    print('synthetic conditioned block proposal PASS',flush=True)
if __name__=='__main__':selftest()
