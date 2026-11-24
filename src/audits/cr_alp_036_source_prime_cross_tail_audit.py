#!/usr/bin/env python3
"""CR-ALP-036 exact Fraction audit. Continuum Gamma/shift identities are human proofs.
No full signed 15,000-mode Schur certificate or Lean/RH conclusion."""
from fractions import Fraction as F
from math import isqrt
a=F(11,10); S=F(177,50); N=15000

def invsqrt_up(n):
    scale=10**12
    v=isqrt(n*scale*scale)+1
    q=F(v,n*scale)
    assert n*q*q>1
    return q

bounds=[]
for k in range(3,13):
    c=1 if k%2 else 2
    lam=F(c,k)+a*a/F(9*k*c)
    assert 0<lam<1
    x=F(11*k,7)/a
    y=F(11*c,7)/a
    A=F(8)*S/(3*a)*(x+lam*y)*(2*a/3)**2
    B=(F(8)*S/(3*a)*(x*x+lam*y*y)+S*(1+lam)/(2*a))*(2*a/3)**3
    bound=A*invsqrt_up(3*N**3)+B*invsqrt_up(5*N**5)
    assert bound<F(1,30000)
    bounds.append(bound)
assert bounds[0]<F(1,100000)
assert sum((q*q for q in bounds),F(0))<F(1,18000)**2
print("PASS: original seven-prime Q_15000 operator on span(w_3,...,w_12) < 1/18000")

for k in range(3,13):
    c=1 if k%2 else 2
    x,y=F(13*k,7),F(13*c,7)
    lam=x*(y*y+F(1,4))/(y*(x*x+F(1,4)))
    for beta in (F(1,2),F(5,2),F(9,2),F(67,2)):
        direct=x/(beta*beta+x*x)-lam*y/(beta*beta+y*y)
        factored=x*(x*x-y*y)*(beta*beta-F(1,4))/(
            (x*x+F(1,4))*(beta*beta+x*x)*(beta*beta+y*y))
        assert direct==factored
        if beta==F(1,2): assert direct==0
print("PASS: J_(1/2) exact high-sine pole-null cancellation algebra")

for k in (3,7,11,12):
    c=1 if k%2 else 2
    lam=F(c,k)+a*a/F(9*k*c)
    j=15001 if k%2 else 15002
    x,y,z=(F(11*t,7)/a for t in (k,c,j))
    M=30000
    tail=z*(x+lam*y)/(12*a*(M-1)**3)
    assert tail<F(1,100000000)
print("PASS: source-specific fixed (j,k) Gamma kernel remainder for M=30000")

# Source-specific physical boundary anchors: D(u)<25, T(u)<4.
assert a+2<16*a**3/F(6)
assert (F(11,7)/a+F(1,2))**2<4
assert (F(22,7)/a+2)**2<25
assert 2*F(6,5)**2/(2*a)<4
assert 2*F(3,5)**2/(a**3/F(6))<4
M=30000
assert 2*M>200**2
b3=F(1,8*M**3)+F(1,16*M**2)
b52=(F(1,M**2)+F(1,M))/800
for j in (15001,15002):
    z=F(11*j,7)/a
    assert 3*z*(5*b3+2*b52)<F(1,100)
print("PASS: boundary-anchor D<25, T<4 and per-index Gamma tail <1/100")
print("CR-ALP-036 EXACT RATIONAL GATES PASS; whole signed Schur OPEN")
