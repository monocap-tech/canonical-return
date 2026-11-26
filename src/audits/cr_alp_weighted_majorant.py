"""Exact weighted interval majorants; candidate bounds, not physical evidence."""
from fractions import Fraction as F


def majorant(center, radius, weights, reserve=F(1,10**18)):
    n=len(center)
    assert n>0 and len(radius)==len(weights)==n and reserve>0
    C=[[F(x) for x in row] for row in center]
    R=[[F(x) for x in row] for row in radius]
    s=list(map(F,weights))
    assert all(len(row)==n for row in C+R) and min(s)>0
    assert all(C[i][j]==C[j][i] and R[i][j]==R[j][i]>=0
               for i in range(n) for j in range(n))
    d=[sum(R[i][j]*s[j]/s[i] for j in range(n))+reserve for i in range(n)]
    A=[[C[i][j]+(d[i] if i==j else 0) for j in range(n)] for i in range(n)]
    # For |E_ij| <= R_ij, S (diag(d)-E) S is strictly diagonally dominant,
    # where S=diag(s). Congruence gives A-(C+E) positive definite.
    assert all(d[i]*s[i]*s[i]-sum(R[i][j]*s[i]*s[j] for j in range(n))>0
               for i in range(n))
    return A


def candidates(center,radius):
    n=len(center)
    ans={'uniform':majorant(center,radius,[1]*n)}
    for factor in (2,4):
        for i in range(n):
            s=[1]*n;s[i]=factor
            ans[f'axis_{i}_scale_{factor}']=majorant(center,radius,s)
    return ans
