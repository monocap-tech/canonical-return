#!/usr/bin/env python3
"""Independent integer interval LDL verifier for CR-ALP-039 12+infinite Schur.

Reads the complete, separate shifted matrix interval output from the
primary source-specific audit. Re-rounds its 10^-25 entries to an
independent 10^-20 integer grid, checks symmetry, and performs fresh
outward-rounded LDL. Rechecks the anchor budgets with exact Fraction.
Does not rerun the analytic derivations or imply Lean/RH certification.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json
ROOT=Path(__file__).resolve().parent
B=10**20
SRC=10**25
assert B<SRC


def ceil_div(a,b):
    assert b>0
    return -((-a)//b)

class IV:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=int(lo)
        self.hi=int(lo) if hi is None else int(hi)
        assert self.lo<=self.hi
    def __add__(self,other):
        return IV(self.lo+other.lo,self.hi+other.hi)
    def __neg__(self):return IV(-self.hi,-self.lo)
    def __sub__(self,other):return self+(-other)
    def __mul__(self,other):
        p=[i*j for i in [self.lo,self.hi] for j in [other.lo,other.hi]]
        return IV(min(p)//B,ceil_div(max(p),B))
    def __truediv__(self,other):
        assert other.lo>0
        vals=[(i*B,j) for i in [self.lo,self.hi] for j in [other.lo,other.hi]]
        return IV(min(a//b for a,b in vals),max(ceil_div(a,b) for a,b in vals))

matrix=ROOT/'CR-ALP-039_full_twelve_shifted_H_interval_matrix.tsv'
rows=matrix.read_text().splitlines()
assert len(rows)==72
M={}
for line in rows:
    i,j,lo,hi=map(int,line.split())
    # Round source interval strictly outward onto an independent grid.
    val=IV(lo//(SRC//B),ceil_div(hi,SRC//B))
    assert (i,j) not in M
    M[i,j]=val


def check(ids):
    for i in ids:
        for j in ids:
            a,b=M[i,j],M[j,i]
            assert max(a.lo,b.lo)<=min(a.hi,b.hi)
    L={};D=[]
    for col,j in enumerate(ids):
        p=M[j,j]
        for k in range(col):
            p=p-L[j,ids[k]]*L[j,ids[k]]*D[k]
        assert p.lo>0, (ids,j,p.lo,p.hi)
        D.append(p)
        for row in ids[col+1:]:
            v=M[row,j]
            for k in range(col):
                v=v-L[row,ids[k]]*L[j,ids[k]]*D[k]
            L[row,j]=v/p
    return [str(v.lo)+'/10^20' for v in D]

p_odd=check([3,5,7,9,11,1])
p_even=check([4,6,8,10,12,2])
# Independently check the formula and summation bounds for original anchors.
S=F(177,50);A=F(11,10);N=15000
assert F(2,N)+F(2,N*N)<F(1,86*86)
assert F(1,3*N**3)<F(1,3000000**2)
assert F(1,5*N**5)<F(1,10**10)**2
for k,g,B0,B1,gt,pt,jt in [
    (1,F(3,5),F(6,5),F(3,10),F(1,170),F(1,30),F(1,25)),
    (2,F(5,3),F(3,5),F(3,5),F(1,120),F(1,23),F(1,19)),
]:
    GG = g*(F(7,10)*B0/F(86)+
          F(121,225)*(F(12001,4000)*B0+F(11,7)*B1)/F(3000000)+
          F(8,945)*A**3*B0/F(10**10))+F(1,10**6)
    vk=F(k)*F(22,7)/(2*A)
    PP = g*F(7,5)*S*B0/F(122)+ (
          (F(8)*S*vk/(3*A)+g*S*(2*B1+A*B0))*(2*A/3)**2)/F(3000000)+ (
          (F(8)*S*(vk*vk+F(1,4))/(3*A)+g*S*B0/2)*(2*A/3)**3)/F(10**10)
    assert GG<gt and PP<pt and GG+PP<jt
    print('INDEPENDENT PASS full original signed anchor',k,'Gamma',float(GG),'prime',float(PP),'joint',float(GG+PP))
print('INDEPENDENT PASS both complete 6x6 exact outward LDL blocks at 10^-20')
print(' Odd pivots:',p_odd)
print(' Even pivots:',p_even)
print(' Verified shifted interval matrix SHA256:',hashlib.sha256(matrix.read_bytes()).hexdigest())
print('ANALYTICAL premises, untested 13..15000 shell, Lean and RH all remain OPEN')
