#!/usr/bin/env python3
"""Exact rational AE constants and negative Neumann-surrogate witness.
Analytic resolvent comparison and T1D-Q's arithmetic norm are imported,
not proved by this script. A negative surrogate is NOT a negative Weil test.
"""
from fractions import Fraction as F
from math import isqrt, factorial

class I:
    def __init__(self,l,h=None): self.l=F(l); self.h=F(l if h is None else h); assert self.l<=self.h
    def __add__(self,o):
        o=iv(o); return I(self.l+o.l,self.h+o.h)
    __radd__=__add__
    def __neg__(self): return I(-self.h,-self.l)
    def __sub__(self,o): return self+-iv(o)
    def __rsub__(self,o): return iv(o)+-self
    def __mul__(self,o):
        o=iv(o); x=[a*b for a in (self.l,self.h) for b in (o.l,o.h)];return I(min(x),max(x))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=iv(o); assert o.l>0; return self*I(1/o.h,1/o.l)
    def round(self,d=10**12): return I(F(self.l*d//1,d),F(-((-self.h*d)//1),d))
    def __repr__(self): return str((float(self.l),float(self.h)))
def iv(x):return x if isinstance(x,I) else I(x)
def log1(x,N=60):
    assert 1<=x<=4
    z=(x-1)/(x+1); s=sum((2*z**(2*k+1)/F(2*k+1) for k in range(N)),F(0))
    return I(s,s+2*z**(2*N+1)/(F(2*N+1)*(1-z*z)))
log2=log1(F(2))
def logq(x):
    x=F(x); assert x>=1
    j=0
    while x>=2:x/=2;j+=1
    return (j*log2+log1(x)).round()
def atan(x,N):
    s=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(N)),F(0)); e=x**(2*N+1)/F(2*N+1)
    return I(min(s,s+(-1)**N*e),max(s,s+(-1)**N*e))
pi=(16*atan(F(1,5),24)-4*atan(F(1,239),10)).round()
assert pi.l>3
lgpi=I(logq(pi.l).l,logq(pi.h).h)
H=sum((F(1,k) for k in range(1,1025)),F(0))
gamma=(I(H)-10*log2-I(F(1,2050),F(1,2048))).round()
m0=(-gamma-pi/2-3*log2-lgpi).round()
assert m0.l> -F(43,8) and m0.h< -F(537,100)
# m_Gamma(n*pi/4) increasing; use x=n*pi/8>3n/8.
N=32768; x=F(3*N,8); K=32; t=F(K)+F(1,4)
s=sum((x*x/((F(j)+F(1,4))*((F(j)+F(1,4))**2+x*x)) for j in range(K)),F(0))
# Decreasing positive tail >= integral from t to infinity = .5 log(1+(x/t)^2)>log(x/t).
lower=-F(43,8)+s+logq(x/t).l
assert lower>F(203,25)+F(1,8)
# Exact original prime census.
def prime(n):return n>=2 and all(n%j for j in range(2,isqrt(n)+1))
pp={}
for p in range(2,55):
    if prime(p):
        n=p
        while n<=54: pp[n]=p; n*=p
assert sorted(pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]
def sqrtq(n,d=10**12):
    r=isqrt(n*d*d);return I(F(r,d),F(r+1,d))
def trig(x,odd):
    # Taylor at midpoint, ordinary derivative remainder and 1-Lipschitz radius.
    m=(x.l+x.h)/2; rad=(x.h-x.l)/2; assert abs(m)<7
    s=sum(((-1)**k*m**(2*k+odd)/factorial(2*k+odd) for k in range(23)),F(0))
    e=abs(m)**(45+odd)/factorial(45+odd)+rad
    return I(s-e,s+e).round()
# Exact pole-null source in the form completion: f(y)=cos(pi*y/2)-c, y=x+2.
w=pi/2; c=(I(1)/(I(1)+pi*pi)).round(); xx=(w/2)*(w/2)
# Gamma diagonal at frequency w: first 64 exact intervals plus decreasing tail.
sm=I(0)
for j in range(64):
    t=F(j)+F(1,4); sm+=xx/(t*(I(t*t)+xx))
sm=sm.round(); tail=xx.h/(2*(F(64)-F(3,4))**2)
mw=(m0+sm+I(0,tail)).round()
# Shift autocorrelation for f(y) on [0,4], using sin/cos at ws.
A=I(0)
for n,p in pp.items():
    h=logq(n); wt=logq(p)/sqrtq(n); ph=(w*h).round()
    si=trig(ph,1); co=trig(ph,0); L=4-h
    corr=L*co/2-si/(2*w)+2*c*si/w+c*c*L
    A+=(2*wt*corr).round()
qN=2*mw+4*c*c*m0-A
mass=2+4*c*c
ray=qN/mass
assert -F(13300,100000)<ray.l<=ray.h<-F(13291,100000)
print('AE EXACT-RATIONAL AUDIT PASS')
print('Gamma lower at n=32768:',float(lower))
print('Original q high margin >1/8; negative index <=32768 (<=16384 per parity), conditional on Q')
print('Neumann-only surrogate pole-null Rayleigh interval:',ray,'< -13/100')
print('Negative surrogate is not a negative original Weil direction.')
# Recover original Gamma via the exact positive boundary corrections.
# j=0 vanishes by the pole constraint. For j>=1, exp(-4a)<=exp(-10)<1/10000.
assert sum((F(10)**j/factorial(j) for j in range(20)),F(0))>10000
corr=I(0)
for j in range(1,64):
    a=F(2*j)+F(1,2); v=I(1)/(I(a*a)+w*w)-c/(a*a)
    assert v.l>0
    corr+=(2*a*a*I(F(9999,10000),1)*v*v).round()
tail=F(1,2)/(F(64)-F(3,4))
full=(qN+corr+I(0,tail))/mass
assert F(2534068,100000000)<full.l<=full.h<F(2931959,100000000)
assert full.l>F(1,40) and full.h<F(3,100)
print('Original signed Weil form on same pole-null trial: Rayleigh interval:',full,'in (1/40,3/100)')
