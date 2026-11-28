#!/usr/bin/env python3
"""CR-ALP-035 exact arithmetic audit of high-sine pole projection.

This checks source-specific rational upper/lower enclosures for the two
pole-tail Grams at N=15000, and finite-dimensional exact algebra for the
rank-two projector and nested shell decomposition. It does NOT prove
continuum Parseval/regularity or the signed Weil Schur inequality, and does
not evaluate the original 15,000-dimensional form or certify RH.
"""
from fractions import Fraction as F
from math import factorial

A = F(11,10)
x = A/2
N = 15000

# Taylor-plus-geometric tails: 0<x<1. For cosh, after the x^4/4! term the next ratio is x^2/(5*6);
# for sinh, after the x^5/5! term the next ratio is x^2/(6*7).
cosh_low = 1+x*x/2
cosh_up = cosh_low + (x**4/F(factorial(4))) / (1-x*x/F(5*6))
sinh_low = x
sinh_up = x+x**3/F(factorial(3)) + (x**5/F(factorial(5))) / (1-x*x/F(6*7))
assert cosh_low < cosh_up < F(6,5)
assert sinh_low < sinh_up < F(3,5)

# Elementary integral test on the parity-resolved reciprocal-square tails:
# 1/(2J) <= sum_{k>=0} 1/(J+2k)^2 <= 1/(2J)+1/J^2.
# For nu_j>=1, (16/25)/nu_j^2 <= nu_j^2/(nu_j^2+1/4)^2 <=1/nu_j^2.
# With 3<pi<22/7, the genuine pole-tail sine Gram scalars are enclosed.

def physical_pole_tail_bracket(J, amplitude_lo, amplitude_hi):
    lo = (F(4)/A) * amplitude_lo**2 * F(16,25) * \
         (2*A/F(22,7))**2 * F(1,2*J)
    hi = (F(4)/A) * amplitude_hi**2 * (2*A/F(3))**2 * \
         (F(1,2*J)+F(1,J*J))
    assert lo>0 and hi>lo
    return lo,hi

G_even_lo,G_even_hi = physical_pole_tail_bracket(N+1,cosh_low,F(6,5))
G_odd_lo,G_odd_hi = physical_pole_tail_bracket(N+2,sinh_low,F(3,5))
assert F(1,20000)<G_even_lo<G_even_hi<F(1,10000)
assert F(1,90000)<G_odd_lo<G_odd_hi<F(1,40000)
print('PASS pole-tail Grams at N=15000:')
print('  even (cosh):',float(G_even_lo),'< G <',float(G_even_hi))
print('  odd  (sinh):',float(G_odd_lo),'< G <',float(G_odd_hi))
print('  reciprocal upper bounds:',float(1/G_even_lo),float(1/G_odd_lo))

# Finite-dimensional EXACT regression (not the zeta operator!).
# The source pole vectors have opposite parity, and each retains high-frequency
# components after finite sine cutoff. This finite toy isolates the algebra.
L=14
low=4
b_even=[F(1,j) if j%2==1 else F(0) for j in range(1,L+1)]
b_odd =[F(1,j) if j%2==0 else F(0) for j in range(1,L+1)]
h=[F((7*j*j+3*j+5)%17-8,23+j) for j in range(1,L+1)]

def dot(u,v):
    return sum((a*b for a,b in zip(u,v)),F(0))

def norm2(v): return dot(v,v)

def tail(v,n): return [F(0) if j<n else v[j] for j in range(len(v))]

def project(n,v):
    Rv=tail(v,n)
    b1,b2=tail(b_even,n),tail(b_odd,n)
    g1,g2=norm2(b1),norm2(b2)
    assert g1>0 and g2>0
    return [Rv[j]-dot(b1,Rv)/g1*b1[j]-dot(b2,Rv)/g2*b2[j]
            for j in range(L)]

Qlow=project(low,h)
assert Qlow[:low]==[F(0)]*low
assert dot(Qlow,b_even)==0 and dot(Qlow,b_odd)==0
assert project(low,Qlow)==Qlow
# Self-adjointness check for independent test vector.
g=[F((5*j*j+2)%13-6,19+j) for j in range(1,L+1)]
assert dot(project(low,h),g)==dot(h,project(low,g))
# Exact Pythagorean pole correction (orthogonal because opposite parity).
Rh=tail(h,low)
b1,b2=tail(b_even,low),tail(b_odd,low)
g1,g2=norm2(b1),norm2(b2)
d1,d2=dot(Rh,b1),dot(Rh,b2)
assert norm2(Qlow)==norm2(Rh)-d1*d1/g1-d2*d2/g2
# Q_M-Q_N is a nested shell projection; exact positivity/idempotence.
M=4
NN=9
shell=[a-b for a,b in zip(project(M,h),project(NN,h))]
assert [a-b for a,b in zip(project(M,shell),project(NN,shell))]==shell
assert dot(shell,project(NN,h))==0
assert norm2(project(M,h))==norm2(shell)+norm2(project(NN,h))
assert norm2(shell)>=0
print('PASS exact rank-two high-sine orthogonal projection, parity correction,')
print('     residual variance identity, and nested-shell Pythagoras')
print('CR-ALP-035 EXACT ARITHMETIC/PROJECTION REGRESSION PASS')
print('Continuum projection, original Weil operator coupling and RH not certified.')
