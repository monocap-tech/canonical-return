#!/usr/bin/env python3
"""CR-INT-T1C pure rational, OUTWARD analytical interval audit.

Replicates the elementary estimates for an explicit C-infinity compact bump:
 b_h(x)=(1-|x|/h)_+, rho_e(x)=Z_e^{-1} exp(-1/(1-(x/e)^2))
 for |x|<e (and 0 outside), phi=b_h*rho_e.
 h=49*ell/100, e=ell/100, ell=1/10000.

The enclosures are analytic inequalities, not quadrature or Lean proofs.
The fixed four-cell source uses translated copies phi(x-j log3).

All interval endpoints are exact rational numbers. Logs and pi are
bounded by geometric/alternating series with explicit remainder bounds;
Euler's constant uses elementary harmonic-number bounds.
"""
from fractions import Fraction as F
from math import isqrt

class I:
    def __init__(self, lo, hi=None):
        self.lo=F(lo)
        self.hi=F(hi) if hi is not None else F(lo)
        assert self.lo<=self.hi
    @staticmethod
    def cv(x):return x if isinstance(x,I) else I(x)
    def __add__(self,b):
        b=self.cv(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+(-self.cv(b))
    def __rsub__(self,b):return self.cv(b)+(-self)
    def __mul__(self,b):
        b=self.cv(b)
        a=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=self.cv(b)
        assert b.lo*b.hi>0
        return self*I(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return self.cv(b)/self
    def __pow__(self,n):
        assert n>=0 and isinstance(n,int)
        res=I(1)
        for _ in range(n):res*=self
        return res
    def __repr__(self):return f'[{float(self.lo):.15g}, {float(self.hi):.15g}]'

def log_rat(x,N=50):
    """Exact rational atanh series enclosure for log(x), x>=1."""
    x=F(x); assert x>=1
    z=(x-1)/(x+1)
    if z==0:return I(0)
    v=z;ss=F(0)
    for k in range(N):
        ss+=2*v/F(2*k+1)
        v*=z*z
    remainder=2*v/(F(2*N+1)*(1-z*z))
    return I(ss,ss+remainder)

def atan_pos(x,N=23):
    """Alternating, exact rational arctan bounds 0<x<1."""
    x=F(x);assert 0<x<1
    ss=F(0)
    for k in range(N):ss+=(-1)**k*x**(2*k+1)/F(2*k+1)
    rem=x**(2*N+1)/F(2*N+1)
    return I(ss,ss+rem) if N%2==0 else I(ss-rem,ss)

def sqrt_int(n,digits=30):
    D=10**digits
    lower=isqrt(n*D*D)
    return I(F(lower,D), F(lower+1,D))

def out(x,places=6):
    """Outward finite decimal strings (integer-only; rigorously enclosing)."""
    scale=10**places
    def floor_div(a,b):return a//b
    lo=floor_div(x.lo.numerator*scale,x.lo.denominator)
    hi=-floor_div(-x.hi.numerator*scale,x.hi.denominator)
    def fmt(z):
        sg = '-' if z<0 else ''
        q,rem=divmod(abs(z),scale)
        return sg+str(q)+('.'+str(rem).zfill(places) if places else '')
    return '['+fmt(lo)+', '+fmt(hi)+']'

ELL=F(1,10000)
h=F(49,100)*ELL
EPS=F(1,100)*ELL
j_lo=F(2,3)*h*(1-EPS/h)**3
j_hi=F(2,3)*h
# Norm derivative by Young: ||phi'||^2 <= ||b_h'||^2=2/h.
q_ell_sq=(F(2)/h)/j_lo*ELL**2
assert q_ell_sq==F(30625,2304)

sqrt3=sqrt_int(3)
log2=log_rat(2,40)
log54=log_rat(F(5,4),35)
log10=3*log2+log54
log3=log_rat(3,85)
log20000=13*log2+4*log54
pi=16*atan_pos(F(1,5),29)-4*atan_pos(F(1,239),17)
log4pi=4*log2-I(log_rat(4/pi.hi,55).lo,log_rat(4/pi.lo,55).hi)
# Harmonic gamma enclosure: 1/(2(n+1)) < H_n-log(n)-gamma < 1/(2n).
n=1000
H_n=sum((F(1,k) for k in range(1,n+1)),F(0))
logn=3*log10
gamma=I(H_n,H_n)-logn-I(F(1,2*(n+1)),F(1,2*n))
const=log4pi+gamma
# x coth(x) lies in [1,1+x^2] for x=ell/2.
logcoth=log20000+I(0,ELL*ELL/4)
g0_over_j=logcoth-const+I(-ELL,q_ell_sq/2)

def invpow3half(k):
    if k%2==0:return I(F(1,3**(k//2)))
    return I(1)/((3**(k//2))*sqrt3)

def w_k_lower(k):
    # w(k*log3+ell) >= 3^(-k/2)*(1-ell/2) / (1-3^(-2k)*(1-2ell))
    return invpow3half(k)*I(1-ELL/2)/(I(1)-F(1,3**(2*k))*I(1-2*ELL))

def w_k_upper(k):
    # w(k*log3-ell) <= 3^(-k/2)/(1-ell/2) / (1-3^(-2k)/(1-2ell))
    return invpow3half(k)/(I(1-ELL/2))/(I(1)-F(1,3**(2*k))/I(1-2*ELL))

# Each A_k/J is >= (3h/2)*w(kd+ell), <= ell*w(kd-ell).
a=[]
for k in range(1,4):
    low=(F(3,2)*h*w_k_lower(k)).lo
    high=(ELL*w_k_upper(k)).hi
    a.append(I(low,high))

# Generalized eigenvalues equal g0/J + prime + archimedean, simplification
# from the verified T1B exact X^* T X pencil.
a1,a2,a3=a
mu_plus=(g0_over_j+(41-5*sqrt3)*log3/73
         +(-83+76*sqrt3)*a1/146
         +(15+32*sqrt3)*a2/73
         +(-33-12*sqrt3)*a3/146)
mu_minus=(g0_over_j+(41+5*sqrt3)*log3/73
          +(83+76*sqrt3)*a1/146
          +(15-32*sqrt3)*a2/73
          +(33-12*sqrt3)*a3/146)
# The difference should cancel the common g0 (use direct expression to avoid wide dependence).
split_naive=(10*sqrt3*log3+(83*a1-64*sqrt3*a2+33*a3))/73
# Correlated sharper lower: w(t+log3)<w(t)/sqrt3, hence A2<A1/sqrt3.
split_reinforced_lower=((10*sqrt3*log3+19*a1+33*a3)/73).lo
split=I(max(split_naive.lo, split_reinforced_lower),split_naive.hi)
# adjacent min phase gap = (11-4sqrt3)/22 times splitting
adjmin=(11-4*sqrt3)/22*split

print('ANALYTICAL RATIONAL INTERVAL AUDIT — exact fractions (no quadrature)')
print('ell=',ELL, ' h=',h,' eps=',EPS)
print('J / ell in',out(I(j_lo/ELL,j_hi/ELL),8))
print('J in',out(I(j_lo,j_hi),12))
print('ell^2 * ||phi\' ||^2/J <= ',q_ell_sq)
for nm,v in [('pi',pi),('log3',log3),('Euler gamma',gamma),('log(4pi)+gamma',const),('log coth(ell/2)',logcoth),('g0/J',g0_over_j)]:
 print(nm,out(v,8))
print('g0 in',out(g0_over_j*I(j_lo,j_hi),12))
for k,v in enumerate(a,1):
    print('A%d/J'%k,out(v,10), ' ; absolute A%d'%k,out(v*I(j_lo,j_hi),13))
for nm,v in [('mu_plus',mu_plus),('mu_minus',mu_minus),('mu_minus-minus-mu_plus',split),('phase minimum adjacent gap',adjmin)]:
 print(nm,out(v,8))
assert mu_plus.lo>7
assert mu_minus.lo>7
assert split.lo>F(1,4)
assert adjmin.lo>F(1,25)
print('PASS: positive restricted Weil form, \u0394selected >1/4, and adjacent phase gaps >1/25; NO global positivity implication.')