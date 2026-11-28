#!/usr/bin/env python3
"""CR-INT-T1D-R: exact-rational archimedean floor and compact Schur-gate audit.

Certifies arithmetic bounds for m_Gamma(0), the 24-prime-power
finite-window zero-frequency symbol, the conditional operator constant,
and a finite algebra example for the Galerkin tail lemma.

The mathematical facts (digamma series monotonicity, compact log-Sobolev
embedding, representation theorem, compact spectral reduction) are *analytic*
inputs of the research note, not proved by this script. No low-mode spectral
positivity, operator identity with semilocal Sonin, Lean, or RH proof.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp


class J:
    def __init__(self,lo,hi=None):
        self.lo=F(lo); self.hi=F(hi) if hi is not None else F(lo)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=o if isinstance(o,J) else J(o)
        return J(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self): return J(-self.hi,-self.lo)
    def __sub__(self,o): return self+(-o if isinstance(o,J) else J(-o))
    def __mul__(self,o):
        o=o if isinstance(o,J) else J(o)
        q=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return J(min(q),max(q))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=o if isinstance(o,J) else J(o)
        assert o.lo>0
        return self*J(1/o.hi,1/o.lo)
    def show(self,n=12):
        d=10**n
        lo=self.lo.numerator*d//self.lo.denominator
        hi=-((-self.hi.numerator*d)//self.hi.denominator)
        def p(k):
            sign='-' if k<0 else '';k=abs(k)
            return sign+str(k//d)+'.'+str(k%d).zfill(n)
        return '['+p(lo)+', '+p(hi)+']'


def log_q(q,N=90):
    q=F(q)
    assert q>=1
    z=(q-1)/(q+1);t=z;tot=F(0)
    for k in range(N):
        tot+=2*t/F(2*k+1); t*=z*z
    return J(tot,tot+2*t/(F(2*N+1)*(1-z*z)))


def atan_q(x,N):
    x=F(x);assert 0<x<1
    z=F(0)
    for k in range(N):z+=((-1)**k)*x**(2*k+1)/F(2*k+1)
    e=x**(2*N+1)/F(2*N+1)
    return J(z,z+e) if N%2==0 else J(z-e,z)


def root_q(n,d=10**18):
    m=isqrt(n*d*d)
    return J(F(m,d),F(m+1,d))


def prime_power(n):
    for p in range(2,n+1):
        if n%p:continue
        if any(p%r==0 for r in range(2,isqrt(p)+1)):continue
        m=n
        while m%p==0:m//=p
        if m==1:return p
    return None


pi=16*atan_q(F(1,5),30)-4*atan_q(F(1,239),18)
assert F(31415,10000)<pi.lo<pi.hi<F(31418,10000)
log2=log_q(2)
log_pi=J(log_q(pi.lo).lo,log_q(pi.hi).hi)
H=sum((F(1,k) for k in range(1,1025)),F(0))
gamma=J(H)-10*log2-J(F(1,2*1025),F(1,2*1024))
assert F(5772,10000)<gamma.lo<gamma.hi<F(5773,10000)
assert pi.hi/2<F(15709,10000) and (3*log2).hi<F(208,100)
assert log_pi.hi<F(1146,1000)

# Exact published digamma identity:
# psi(1/4) = -gamma - pi/2 - 3 log2;
# mGamma(0)=psi(1/4)-log pi, and is the *global minimum* by
# Re psi(a+iy)-psi(a)=sum_n y^2/[ (n+a)((n+a)^2+y^2) ] >=0.
mg0=-gamma-pi/2-3*log2-log_pi
assert mg0.lo>-F(43,8)
assert mg0.hi<-5
assert F(6)+mg0.lo>F(5,8)
pp={n:prime_power(n) for n in range(2,55) if prime_power(n)}
assert len(pp)==24 and list(pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]

def log_int(p):
    j=p.bit_length()-1
    return j*log2+log_q(F(p,2**j))

S=J(0)
for n,p in pp.items():S+=log_int(p)/root_q(n)
assert F(1219,100)<S.lo<S.hi<F(122,10)
fullzero=mg0-2*S
assert -30<fullzero.lo<fullzero.hi<-29
assert F(43,8)+F(203,25)==F(2699,200)

# The negative zero-frequency value of the formal truncated whole-line
# multiplier does NOT give a negative pole-null compact source test.
# T1D-Q's |A|<203/25 and T1D-I's |R|<155/182 imply |C|<15,
# conditional on that separately derived T1D-I residual enclosure.
A_upper=F(203,25)
R_upper=F(155,182)
assert 6+A_upper+R_upper<15

# A finite symbolic sanity check of the Galerkin tail inequality.
D=sp.diag(1,4,16,64)
C=sp.Matrix([[2,1,-1,0],[1,0,1,1],[-1,1,2,-2],[0,1,-2,3]])
E=sp.diag(1,sp.Rational(1,2),sp.Rational(1,4),sp.Rational(1,8))
B=E*C*E
P=sp.diag(1,1,0,0)
Diff=B-P*B*P
rowbound=max([sum(abs(C[i,j]) for j in range(4)) for i in range(4)])
assert rowbound==6
frob_sq=sum(Diff[i,j]**2 for i in range(4) for j in range(4))
epsilon=2*F(rowbound,4)+F(rowbound,16)
assert frob_sq<epsilon*epsilon

print('PASS T1D-R: global Gamma multiplier minimum at xi=0 (monotonicity: analytic input)')
print('CERTIFIED mGamma(0)',mg0.show(),'> -43/8; hence Gamma+6 >= 5/8')
print('CERTIFIED genuine prime sum S',S.show(),'(24 prime powers)')
print('CERTIFIED formal whole-line truncated symbol at xi=0:',fullzero.show(),'in (-30,-29)')
print('CERTIFIED whole pole-null source semibound: q(y)> -(2699/200)||y||^2 (analytic Gamma floor + Q norm)')
print('PASS: conditional ||C_F||<15 from Q prime bound + I cross penalty')
print('PASS: Galerkin tail norm inequality, finite symbolic sanity check')
print('ANALYTIC: compact Birman-Schwinger criterion, smooth form core, parity decomposition')
print('OPEN: top eigenvalue for either parity, all-F Schur positivity, full pole directions, RH')
