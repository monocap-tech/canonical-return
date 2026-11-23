#!/usr/bin/env python3
"""CR-INT-T1D-AA: exact-rational Fourier-localized second Schur certificate.

Checks the 24-prime-power whole-line multiplier bound, source H1
budgets, concentration leakage, low/high mixed functional inequalities,
parity-specific rank-two/second-rank-four penalty bounds, and the
strict conditioned second-Schur lower margin at B=2^21, theta=2^-38.

Analytic inputs, not proved by this script: T1D-F Gamma symbol
|m_Gamma| <= log(1+|xi|)+40; T1D-M/O/Q Fourier concentration and
Gamma lower inequality; T1C selected derivative control; T1D-W/X/Y
source couplings, positive 2x2 Schur matrices, and source/core facts.
No low-frequency complement positivity, complete spectral top, Lean or RH.
"""
from fractions import Fraction as F
from math import isqrt, factorial

class Iv:
    def __init__(self, lo, hi=None):
        self.lo=F(lo); self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=o if isinstance(o,Iv) else Iv(o)
        return Iv(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __mul__(self,o):
        o=o if isinstance(o,Iv) else Iv(o)
        a=[x*y for x in (self.lo,self.hi) for y in (o.lo,o.hi)]
        return Iv(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=o if isinstance(o,Iv) else Iv(o)
        assert o.lo>0
        return self*Iv(F(1,o.hi),F(1,o.lo))
    def __str__(self):
        return f'({float(self.lo):.12f}, {float(self.hi):.12f}) [diagnostic only]'

def logr(r,N=80):
    r=F(r); assert r>=1
    a=(r-1)/(r+1); term=a; total=F(0)
    for j in range(N):
        total+=2*term/F(2*j+1)
        term*=a*a
    return Iv(total,total+2*term/(F(2*N+1)*(1-a*a)))

def root(n,scale=10**15):
    m=isqrt(n*scale*scale)
    assert m*m<=n*scale*scale<(m+1)*(m+1)
    return Iv(F(m,scale),F(m+1,scale))

def prime_base(n):
    for p in range(2,n+1):
        if n%p: continue
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        x=n
        while x%p==0:x//=p
        if x==1:return p
    return None

pp={n:prime_base(n) for n in range(2,55) if prime_base(n)}
assert list(pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]
log2=logr(2)
def logint(p):
    k=p.bit_length()-1
    return Iv(k)*log2+logr(F(p,2**k))
S=sum((logint(p)/root(n) for n,p in pp.items()),Iv(0))
assert 2*S.hi<25
assert 2*S.hi>24
print('PASS: 24 original prime powers; whole-line prime symbol |a(xi)|<=2S<25; 2S:',2*S)

# New work uses whole-line mixed source representer K v=(Gamma-A_R)v,
# rather than zero-extending the physically projected representer.
# T1D-F |m_Gamma|<=log(1+|xi|)+40; |a_R|<25.
# Therefore |k(xi)| <=log(1+|xi|)+65.
# Pointwise log(1+x)<=x gives ||K v|| <=(65+||v'||/||v||)||v||.
# For x>=B, (log(1+x)+65)/x is decreasing; if log(1+B)<15,
# then ||(K v)hat_{>B}||<=80/B *||v'||.
B=2**21
# Exact positive Taylor lower bound: exp(15)> 1+2**21.
exp15_lower=sum((F(15)**k/F(factorial(k)) for k in range(41)),F(0))
assert exp15_lower>1+B
assert log2.lo>F(693,1000)
print('PASS: exp(15)>2^21+1; uniform high-frequency multiplier envelope <=80/B')

# T1C: ell^2*||phi'||^2/J <=30625/2304, ell=1/10000.
# The normalized selected combinations use pairwise disjoint bumps,
# so the same normalized derivative budget applies to u_+ and u_-.
ell=F(1,10000)
r_u=37000
assert F(30625,2304)/ell**2<r_u*r_u
# T1D-V/W/X/Y: four normalized central/satellite derivative ratios
# <=38 (even), <=46 (odd), and q(u,v) bounds <=27/100.
r_v=46
Btheta=F(1,2**38)
root_theta=F(1,2**19)
assert root_theta**2==Btheta

# Cauchy-Schwarz frequency split, for w with low-band mass <=theta:
# |q(v,w)|/||w|| <= (65+r_v)*sqrt(theta)+80*r_v/B.
# Selected |q(u_sigma,w)|/||w|| <=(65+r_u)*sqrt(theta)+80*r_u/B.
U=F(65+r_u,2**19)+F(80*r_u,B)
V=F(65+r_v,2**19)+F(80*r_v,B)
assert U<F(3,2),(U,float(U))
assert V<F(1,500),(V,float(V))
print('PASS: exact whole-line Fourier split selected cross U=',U,'(~%.12f)'%float(U),'<3/2')
print('PASS: exact whole-line Fourier split source cross V=',V,'(~%.12f)'%float(V),'<1/500')

# Original selected mu_sigma>7, T1D-W/Y selected-to-new-source:
# max over four unit vectors |q(u_sigma,v_i)|<27/100.
# Each b_i(w)=q(v_i,w)-q(v_i,u_sigma)q(u_sigma,w)/mu_sigma.
c=F(27,100)
bound_b=F(1,500)+(c/F(7))*F(3,2)
assert bound_b<F(3,50)
# Parity: only ONE selected direction acts per sector, giving R[w]
# <=(3/2)^2/7 ||w||^2. The 2x2 Schur matrices in L have
# min eigenvalues >1/70 (even), >1/18 (odd) from T1D-Z.
R_bound=F(3,2)**2/F(7)
T_even=F(70)*2*F(3,50)**2
T_odd=F(18)*2*F(3,50)**2
assert R_bound==F(9,28) and T_even==F(63,125)
assert T_odd<T_even
print('PASS: |b_i(w)|<3/50; parity R[w]<9/28; second Schur T_even[w]<63/125, odd smaller')

# T1D-Q/M/O: on the original smooth prolate-conditioned F source,
# gamma >=(1-theta)(log B-49/10)-6*theta. Combined q>=gamma-203/25.
weillower=(1-Btheta)*(21*F(693,1000)-F(49,10))-6*Btheta-F(203,25)
assert weillower>F(153,100)
finalmargin=F(153,100)-R_bound-T_even
assert finalmargin==F(1233,1750)>F(7,10)
print('PASS: genuine full original q(w)>153/100 ||w||^2 on B=2^21, theta=2^-38 conditioned source')
print('CERTIFIED: q(w)-R[w]-T[w] > 1233/1750 ||w||^2 > 7/10 ||w||^2')
print('No four additional mixed-null conditions needed; selected-to-remainder and second Schur controlled together.')
print('OPEN: complement excluded by prolate projection, entire low-frequency remainder, global Weil positivity, RH.')
