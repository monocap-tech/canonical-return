#!/usr/bin/env python3
"""CR-INT-T1D-H independent exact-rational analytic audit.

Certifies the Fourier-log variance of the actual mollified T1C smooth profile
via an exact triangular reference Mellin transform and a compact-mollifier
Fourier perturbation inequality. Certifies nonzero finite-width mixed Weil
couplings with one explicit, pole-null, orthogonal translated witness.

No numerical quadrature, no false certification of the whole residual norm,
no infinite-dimensional complement positivity, Lean, or RH.
"""
from fractions import Fraction as F
from math import isqrt

class I:
    def __init__(self,lo,hi=None):
        self.lo=F(lo)
        self.hi=F(hi) if hi is not None else F(lo)
        assert self.lo<=self.hi
    @staticmethod
    def c(x): return x if isinstance(x,I) else I(x)
    def __add__(self,o):
        o=self.c(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-self.c(o))
    def __rsub__(self,o):return self.c(o)-self
    def __mul__(self,o):
        o=self.c(o); z=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return I(min(z),max(z))
    __rmul__=__mul__
    def inv(self):
        assert self.lo*self.hi>0
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,o):return self*self.c(o).inv()
    def __rtruediv__(self,o):return self.c(o)*self.inv()
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        if n==0:return I(1)
        if n==1:return self
        if n%2==0 and self.lo<0<self.hi:
            return I(0,max(abs(self.lo),abs(self.hi))**n)
        if n%2==0:return I(min(self.lo**n,self.hi**n),max(self.lo**n,self.hi**n))
        return I(self.lo**n,self.hi**n)

def sqrt_int(n, digits=30):
    D=10**digits; z=isqrt(n*D*D)
    assert z*z<=n*D*D<(z+1)*(z+1)
    return I(F(z,D),F(z+1,D))

def log_rat(x,N=80):
    x=F(x);assert x>=1
    z=(x-1)/(x+1); t=z;sm=F(0)
    for j in range(N):
        sm+=2*t/F(2*j+1);t*=z*z
    return I(sm,sm+2*t/(F(2*N+1)*(1-z*z)))

def atan_pos(x,N=28):
    x=F(x);assert 0<x<1
    a=sum((F((-1)**k)*x**(2*k+1)/F(2*k+1) for k in range(N)),F(0))
    t=x**(2*N+1)/F(2*N+1)
    return I(a,a+t) if N%2==0 else I(a-t,a)

def outward(z,n=8):
    m=10**n
    lo=z.lo.numerator*m//z.lo.denominator
    hi=-((-z.hi.numerator*m)//z.hi.denominator)
    def s(k):
        neg='-' if k<0 else '';k=abs(k)
        return neg+str(k//m)+'.'+str(k%m).zfill(n)
    return '['+s(lo)+', '+s(hi)+']'

# Exact rational constants, with analytic remainder bounds.
pi=16*atan_pos(F(1,5),29)-4*atan_pos(F(1,239),17)
log2=log_rat(2,90)
log3=log_rat(3,100)
log3over2=log_rat(F(3,2),90)
n=1000
H=sum((F(1,k) for k in range(1,n+1)), F(0))
log10=3*log2+log_rat(F(5,4),90)
gamma=I(H)-3*log10-I(F(1,2*(n+1)),F(1,2*n))
assert F(3)<pi.lo<pi.hi<F(22,7)
assert F(69,100)<log2.lo<log2.hi<F(7,10)
assert F(1,2)<gamma.lo<gamma.hi<F(3,5)

# The triangular profile b_h(x)=(1-|x|/h)_+, h=49/100,
# normalized spectral p(x)=3/(2*pi)*(sin(x)/x)^4, where x=eta*h/2.
# Mellin transform M(z)=3*2^(2-z)*(2^(1-z)-1)/
#                     (2*Gamma(4-z)*cos(pi*z/2)).
# Derivatives at z=0 and z=2 yield exact reference moments.
ref_mean=I(F(11,6))-gamma-3*log2
Vtri=pi*pi/12+I(F(49,36))-2*log2*log2
ref_second=Vtri+ref_mean*ref_mean
weighted_log_square=I(F(3,4))*(1+pi*pi/12-2*log2*log2+(1-gamma)*(1-gamma))
assert Vtri.lo>F(6,5)
assert weighted_log_square.lo>0

# rho_eps nonnegative symmetric probability supported |t|<=eps,
# hence 0 <= 1-|rhohat(eta)|^2 <= eps^2 eta^2.
# Actual spectral law p_smooth proportional to p_triangle*|rhohat|^2.
# eps/h=1/49, x=eta*h/2. Exact triangular E[x^2]=3/4.
r=F(1,49)
D0=F(3,1)*r*r        # E[1-w] <= 3/2401
D1=4*r*r*(weighted_log_square+I(F(3,4)))/2
D2=4*r*r*weighted_log_square
mean=(ref_mean+I(-D1.hi,D1.hi))/I(1-D0,1)
second=I((ref_second-D2).lo,(ref_second/I(1-D0)).hi)
V=second-mean**2
assert V.lo>F(6,5) and V.hi<F(5,4)
print('CERTIFIED TRIANGULAR LOG VARIANCE',outward(Vtri,10))
print('CERTIFIED ACTUAL SMOOTH V_v',outward(V,8))
print('CERTIFIED smooth-variance containment: 6/5 < V_v < 5/4')

# Frozen-width nonzero mixed coupling: y=T_log2 e0 / sqrt(J Ny).
# The 3 shifted cells lie inside (-2,2), both Mellin moments vanish,
# and they are disjoint from all 4 cells of u_plus,u_minus.
# At ell=1/10000 the ONLY prime-power coincidence is prime 2,
# by rational ratio checks against every true p^k <=54.
def primepower(z):
    for p in range(2,z+1):
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        a=z
        while a%p==0:a//=p
        if a==1:return True
    return False
primepowers=[z for z in range(2,55) if primepower(z)]
assert len(primepowers)==24
centerdist=[]
for i in range(4):
    for j in range(3):
        k=j-i
        shift=F(2)*F(3)**k
        v=max(shift,1/shift)
        for npp in primepowers:
            rat=max(v/npp,F(npp)/v)
            if (k,npp)==(0,2): assert rat==1
            else: centerdist.append(rat)
assert min(centerdist)==F(27,26),min(centerdist)
ell=F(1,10000)
assert F(1)/(1-ell)<min(centerdist)
# nearest original/shifted center distance is log(3/2).
assert (2*(F(1,5)+F(1,5)**3/3)-ell)>F(2,5)
assert sum(F(4,5)**k/F(__import__('math').factorial(k)) for k in range(5))>2
# Consequently w(t)=exp(-t/2)/(1-exp(-2t))<2 for t>=log(3/2)-ell.

sqrt3=sqrt_int(3);sqrt2=sqrt_int(2)
s=4/sqrt3
Ny=I(F(22,3))
Nplus=(44-16*sqrt3)/3
Nminus=(44+16*sqrt3)/3
kappa2=log2/sqrt2
# All-source absolute Gamma mixed bound <= 2ell ||c||_1||d||_1/sqrt(Nsigma Ny).
# sqrt(Nsig Ny) is enclosed using exact integer-square-root square root of interval
# approximations; sqrt interval by integer isqrt rational endpoints.
def sqrt_frac_interval(z, digits=28):
    D=10**digits
    a=(z.lo.numerator*D*D)//z.lo.denominator
    b=-((-z.hi.numerator*D*D)//z.hi.denominator)
    low=isqrt(a);hi=isqrt(b)+1
    return I(F(low,D),F(hi,D))

coupling={}
for tag,N,c1,C in [
    ('+',Nplus,I(2)+2*(s-1),2-2*s+s*s),
    ('-',Nminus,I(2)+2*(s+1),2+2*s+s*s)
]:
    # cabs sums: plus (1,s-1,s-1,1), minus (1,s+1,s+1,1).
    denom=sqrt_frac_interval(N*Ny)
    prime=kappa2*C/denom
    gamma_error=2*ell*c1*(2+s)/denom
    lower=prime.lo-gamma_error.hi
    upper=prime.hi+gamma_error.hi
    coupling[tag]=(lower,upper,prime,gamma_error)
    print('source',tag,'prime',outward(prime,8),'Gamma magnitude <=',outward(gamma_error,8))
    print('CERTIFIED |q(u_'+tag+',y)| / ||y|| >=',outward(I(lower),8))
    assert lower>(F(21,100) if tag=='+' else F(11,25))

# Because y is an explicit unit vector in F, beta_sigma >= |q(u_sigma,y)|.
# beta^2 are rigorously bounded below at the exact frozen width.
print('CERTIFIED beta_plus^2 > (21/100)^2; beta_minus^2 > (11/25)^2 at ell=1/10000')
# T1C independently certified mu_minus<14.18874262<71/5; combined with
# this rigorous source witness, ||R|| > (11/25)^2/(71/5)=121/8875.
assert F(11,25)**2/F(71,5)==F(121,8875)
print('CONDITIONAL on independently verified T1C mu_minus<71/5: ||R|| > 121/8875.')
print('PASS T1D-H: smooth Fourier-log variance bounded by exact Mellin reference;')
print('explicit smooth pole-null source witness yields nonzero finite-width coupling.')
print('OPEN: actual finite-width beta_+/beta_- numerical enclosures,');print('complement lower bound, completed-domain transfer, Lean, RH.')
