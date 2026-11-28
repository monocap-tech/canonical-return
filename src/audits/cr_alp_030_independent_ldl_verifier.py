#!/usr/bin/env python3
"""Independent Fraction-interval LDL verification of CR-ALP-030 output.

Checks rational enclosure validity, symmetry-compatible interval pairs,
positive physical Gram and strict shifted H-10^-8 M for both parity blocks.
This DOES NOT by itself rederive continuum Gamma/sine Laplace identities.
"""
from fractions import Fraction as F
from pathlib import Path
import json, hashlib

HERE=Path(__file__).resolve().parent
SCALE=10**25
INITIAL=HERE/'CR-ALP-030_initial_rational_intervals_1e25.txt'
HFILE=HERE/'CR-ALP-030_H30000_exact_matrix_intervals.tsv'
assert INITIAL.exists() and HFILE.exists()

class Interval:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        if isinstance(lo,Interval):self.lo,self.hi=lo.lo,lo.hi;return
        self.lo,self.hi=F(lo),F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self, other):
        b=Interval(other);return Interval(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,other):return self+-Interval(other)
    def __rsub__(self,other):return Interval(other)+-self
    def __mul__(self,other):
        b=Interval(other)
        vals=(self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi)
        return Interval(min(vals),max(vals))
    __rmul__=__mul__
    def __truediv__(self,other):
        b=Interval(other)
        assert b.lo>0,"Cannot interval-divide by a nonpositive range"
        return self*Interval(1/b.hi,1/b.lo)

def from_fixed(lo,hi):return Interval(F(int(lo),SCALE),F(int(hi),SCALE))
initial={}
for line in INITIAL.read_text().splitlines():
    name,lo,hi=line.split()
    assert name not in initial
    initial[name]=from_fixed(lo,hi)
assert len(initial)==89

# The analytic kernel is symmetric. Take hulls, not unjustified point values.
entries={}
for line in HFILE.read_text().splitlines():
    i,j,lo,hi=line.split();ij=(int(i),int(j))
    assert ij not in entries
    entries[ij]=from_fixed(lo,hi)
assert len(entries)==32

def symmetric(matrix,i,j):
    a=matrix[(i,j)];b=matrix[(j,i)]
    assert max(a.lo,b.lo)<=min(a.hi,b.hi), (i,j,"incompatible enclosure")
    return Interval(min(a.lo,b.lo),max(a.hi,b.hi))

def H(i,j):return symmetric(entries,i,j)
def M(i,j):return symmetric({(i,j):initial[f'M_{i}_{j}'] for i,j in entries},i,j)

def ldl(order,shift=F(0)):
    n=len(order); piv=[];L={}
    for t,i in enumerate(order):
        v=H(i,i)-shift*M(i,i)
        for k in range(t):v=v-L[(t,k)]*L[(t,k)]*piv[k]
        assert v.lo>0,("Uncertified LDL pivot",order,t,v.lo,v.hi)
        piv.append(v)
        for u in range(t+1,n):
            j=order[u]; x=H(j,i)-shift*M(j,i)
            for k in range(t):x=x-L[(u,k)]*L[(t,k)]*piv[k]
            L[(u,t)]=x/v
    return piv

def mass_ldl(order):
    n=len(order);L={};p=[]
    for t,i in enumerate(order):
        v=M(i,i)
        for k in range(t):v=v-L[(t,k)]*L[(t,k)]*p[k]
        assert v.lo>0,("Physical Gram nonpositive",order,t)
        p.append(v)
        for u in range(t+1,n):
            x=M(order[u],i)
            for k in range(t):x=x-L[(u,k)]*L[(t,k)]*p[k]
            L[(u,t)]=x/v
    return p

orders={'even physical':(3,5,7,1),'odd physical':(4,6,8,2)}
pivots={}
for name,order in orders.items():
    assert all(p.lo>0 for p in mass_ldl(order))
    piv=ldl(order,F(1,10**8))
    pivots[name]=[str(p.lo) for p in piv]
    print(name,'SHIFTED 1e-8 exact Fraction LDL: PASS')
    for t,p in enumerate(piv):
        print(f'  pivot {t+1} > {float(p.lo):.14g}')

report={
  'stage':'CR-ALP-030',
  'status':'EXACT OUTWARD-INTEGER + INDEPENDENT RATIONAL LDL PASS, ANALYTIC CONTINUUM INPUT; NOT LEAN',
  'finite_gamma_terms':30000,
  'basis':'w3,w5,w7,u1 and w4,w6,w8,u2, true physical pole-null Gram',
  'original_signs':'two original poles vanish; all seven paired von Mangoldt prime contributions negative; full original Gamma',
  'physical_lower_bound_strict':'1/100000000',
  'verified_shifted_LDL_pivot_lower_fractions':pivots,
  'input_sha256':hashlib.sha256(INITIAL.read_bytes()).hexdigest(),
  'matrix_intervals_sha256':hashlib.sha256(HFILE.read_bytes()).hexdigest(),
  'limitations':'No 15k finite Schur inequality; no cross-tail residual; no arbitrary support; no original RH; no Lean kernel certification'
}
out=HERE/'CR-ALP-030_certified_eight_mode_result.json'
out.write_text(json.dumps(report,indent=2)+'\n')
print('CR-ALP-030 INDEPENDENT EXACT PHYSICAL 8x8 INTERVAL VERIFIER PASS')
print('Result file:',out)
