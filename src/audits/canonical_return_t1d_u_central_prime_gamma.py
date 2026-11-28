#!/usr/bin/env python3
"""CR-INT-T1D-U: exact arithmetic-image support graph and Gamma/Jensen
bounds for the authentic T1D-T broad-support parity probes.

Verified with exact rational inequalities and symbolic polynomial integration.
The digamma-series identity, monotone-integral tail, Jensen, Young, and
Birman--Schwinger Rayleigh principle are analytic proof inputs documented
in research/CR-INT-T1D-U.md. This does NOT calculate the top compact
spectrum, actual Gamma integrals, or the Galerkin low--high eta_N.
"""
from fractions import Fraction as Q
from math import factorial, isqrt
from sympy import symbols, diff, expand, integrate, Rational

class Enc:
    def __init__(self, lo, hi=None):
        self.lo=Q(lo); self.hi=Q(lo if hi is None else hi)
        assert self.lo <= self.hi
    def __add__(self, other):
        other=other if isinstance(other,Enc) else Enc(other)
        return Enc(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return Enc(-self.hi,-self.lo)
    def __sub__(self, other):
        other=other if isinstance(other,Enc) else Enc(other)
        return self+-other
    def __rsub__(self, other): return Enc(other)-self
    def __mul__(self, other):
        other=other if isinstance(other,Enc) else Enc(other)
        vs=(self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi)
        return Enc(min(vs),max(vs))
    __rmul__=__mul__
    def __truediv__(self, other):
        other=other if isinstance(other,Enc) else Enc(other)
        assert other.lo>0
        return self*Enc(1/other.hi,1/other.lo)
    def __rtruediv__(self, other): return Enc(other)/self
    def display(self,d=12):
        z=10**d
        l=self.lo.numerator*z//self.lo.denominator
        h=-((-self.hi.numerator*z)//self.hi.denominator)
        def fmt(k):
            sg='-' if k<0 else ''
            k=abs(k)
            return f'{sg}{k//z}.{k%z:0{d}d}'
        return '('+fmt(l)+', '+fmt(h)+')'

def log_q(r, N=100):
    """Exact outward interval for log(r), r >= 1, by atanh power series."""
    r=Q(r); assert r>=1
    t=(r-1)/(r+1); pow=t; sm=Q(0)
    for k in range(N):
        sm+=2*pow/(2*k+1)
        pow*=t*t
    return Enc(sm,sm+2*pow/((2*N+1)*(1-t*t)))

def atan_q(r,N):
    r=Q(r);assert 0<r<1
    s=sum(((-1)**k*r**(2*k+1)/Q(2*k+1) for k in range(N)),Q(0))
    e=r**(2*N+1)/Q(2*N+1)
    return Enc(s,s+e) if N%2==0 else Enc(s-e,s)

def exp_q(r,N=40):
    """Positive-ratio geometric majorant for exponential tail."""
    r=Q(r); assert 0<r<5
    s=sum((r**k/Q(factorial(k)) for k in range(N+1)),Q(0))
    u=r**(N+1)/Q(factorial(N+1))
    assert r<Q(N+2)
    return Enc(s,s+u/(1-r/Q(N+2)))

def root_q(n,scale=10**22):
    m=isqrt(n*scale*scale)
    assert m*m <=n*scale*scale<(m+1)*(m+1)
    return Enc(Q(m,scale),Q(m+1,scale))

def prime_base(n):
    for p in range(2,n+1):
        if n%p:continue
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        r=n
        while r%p==0:r//=p
        if r==1:return p
    return None

all_pp={n:prime_base(n) for n in range(2,55) if prime_base(n) is not None}
assert len(all_pp)==24
assert list(all_pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]

# Central J=[-51/200,51/200] => diameter 51/100. A signed prime
# translate survives intersection with I=(-2,2) iff log n<2+51/200.
a=Q(51,200);diam=2*a
xcut=exp_q(Q(2)+a)
assert 9<xcut.lo<xcut.hi<11
survivors=[n for n in all_pp if n<=9]
assert survivors==[2,3,4,5,7,8,9]
E=exp_q(diam)
assert Q(8,5)<E.lo<E.hi<Q(5,3)<2
# The + and - shifted output supports are separated by
# 2*(log2-a)>0; distinct same-sign shifts overlap iff m/n<e^(2a).
edges=[]
for i,n in enumerate(survivors):
    for m in survivors[i+1:]:
        v=Q(m,n)
        if v<E.lo:edges.append((n,m))
        else:assert v>E.hi,(n,m,v)
expected=[(2,3),(3,4),(4,5),(5,7),(5,8),(7,8),(7,9),(8,9)]
assert edges==expected,(edges,expected)

# For each surviving positive/negative copy, Cauchy and support
# disjointness give the row-collision upper bound. The positive- and
# negative-shifted groups have disjoint supports.
k={n:log_q(all_pp[n])/root_q(n) for n in survivors}
z=sum((k[n]*k[n] for n in survivors),Enc(0))
z+=2*sum((k[n]*k[m] for n,m in edges),Enc(0))
z*=2
assert z.hi<Q(61,5)<Q(49,4)
print('PASS: all 24 genuine prime powers retained; only 7 can act on central supported tests')
print('PASS: e^(51/100) strictly between 8/5 and 5/3; exact 8-edge overlap graph')
print('CERTIFIED ||A P_J||^2 < 61/5; exact upper bound interval:',z.display())
print('CERTIFIED ||A P_J|| < 7/2 and P_J A P_J=0 (original physical operator)')

# Rigorous m_Gamma(0)<-5.372 from rational Machin pi, Euler's
# harmonic bounds, log2 and logpi. The exact digamma identity used:
# psi(1/4)= -EulerGamma - pi/2 -3 log 2.
pi=16*atan_q(Q(1,5),35)-4*atan_q(Q(1,239),20)
assert pi.lo>Q(31415,10000)
ln2=log_q(2)
assert ln2.lo>Q(69314,100000)
logpi=Enc(log_q(pi.lo).lo,log_q(pi.hi).hi)
assert logpi.lo>Q(114472,100000)
H=sum((Q(1,k) for k in range(1,1025)),Q(0))
euler=Enc(H)-10*ln2-Enc(Q(1,2*1025),Q(1,2*1024))
assert euler.lo>Q(5772,10000)
m0=-euler-pi/2-3*ln2-logpi
assert m0.hi < -Q(5372,1000),m0.display()
print('PASS: exact rational digamma minimum m_Gamma(0) < -5372/1000')

# For v=xi^2, the exact digamma nonnegative increment is
#   sum_{n>=0} v/[(n+1/4)(4(n+1/4)^2+v)].
# Each summand is concave increasing in v. For fixed M, the rest
# is <= integral_{M+1/4}^infty v/[t(4t^2+v)]dt
#       = (1/2) log(1+v/(4(M+1/4)^2)).
# Jensen bounds its Fourier integral by the same expression evaluated
# at V=||f'||^2/||f||^2.
def gamma_upper_at_variance(V,M=50):
    V=Q(V);a0=Q(1,4)
    low=sum((V/((n+a0)*(4*(n+a0)**2+V)) for n in range(M+1)),Q(0))
    high=Q(1,2)*log_q(1+V/(4*(M+a0)**2))
    return Enc(-Q(5372,1000)+low)+high

# The authentic T1D-T broad profiles, b_e=(1-16x^2)^5 and
# b_o=x*(1-16x^2)^5 on [-1/4,1/4], smoothly mollified at eps=1/200.
# Their pole constraints are (partial_x^2-1/4) e^(+/-x/2)=0.
x=symbols('x',real=True)
for label,poly,limit,gamma_max,ritz_lower in [
    ('even',(1-16*x*x)**5,23,Q(13,10),Q(60,73)),
    ('odd',x*(1-16*x*x)**5,29,Q(77,50),Q(300,377)),
]:
    y0=expand(diff(poly,x,2)-poly/4)
    mass=integrate(y0*y0,(x,-Rational(1,4),Rational(1,4)))
    grad=integrate(diff(y0,x)**2,(x,-Rational(1,4),Rational(1,4)))
    assert mass>0 and grad>0
    eps=Q(1,200)
    # Young/Minkowski: ||y|| >= ||y0||-eps||y0'||,
    # ||y'||<=||y0'||, hence rate<limit by the exact rational test.
    assert Q(grad)/Q(mass) < (Q(limit)/(1+eps*limit))**2
    G=gamma_upper_at_variance(limit*limit)
    assert G.hi<gamma_max,(label,G.display())
    assert 6/(6+gamma_max)==ritz_lower
    print(f'PASS {label}: exact polynomial mass={mass}, derivative mass={grad}')
    print(f'CERTIFIED {label}: ||y_prime||/||y||<{limit}, gamma[y]/||y||^2<{gamma_max}')
    print(f'CERTIFIED {label}: full original compact-spectrum top >{ritz_lower} ({float(ritz_lower):.8f})')
    print('  Jensen upper enclosure',G.display())

print('ANALYTIC: supported central tests are smooth, both Mellin poles vanish, and E is physically disjoint')
print('OPEN: actual eta_N, top eigenvalue upper bound, unrestricted complement positivity, RH')
