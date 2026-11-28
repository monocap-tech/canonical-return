#!/usr/bin/env python3
"""CR-INT-T1D-K independent exact rational / symbolic source audit.

Theorem scope: original signed Weil form on six explicitly constructed,
physically orthonormal, two-Mellin-pole-null smooth directions at ell=1e-4.
All actual analytical Gamma estimates are proved in the companion note;
this script verifies their rational bounds and complete prime masks.
The residual-unmeasured lower bounds additionally import the previous
T1D-I analytic beta^2 inequalities; no Lean, global positivity, or RH.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as S

ELL=F(1,10000)

class Interval:
    def __init__(self,lo,hi=None):
        self.lo=F(lo);self.hi=F(hi) if hi is not None else F(lo)
        assert self.lo<=self.hi
    @staticmethod
    def cv(x):return x if isinstance(x,Interval) else Interval(x)
    def __add__(self,b):
        b=self.cv(b);return Interval(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,b):return self+(-self.cv(b))
    def __rsub__(self,b):return self.cv(b)+(-self)
    def __mul__(self,b):
        b=self.cv(b);ps=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return Interval(min(ps),max(ps))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=self.cv(b);assert b.lo*b.hi>0
        return self*Interval(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return self.cv(b)/self
    def square(self):
        return Interval(0,max(self.lo**2,self.hi**2)) if self.lo<0<self.hi else Interval(min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))

def sqr(n,digits=30):
    q=10**digits;t=isqrt(n*q*q)
    return Interval(F(t,q),F(t+1,q))
def sqrt_interval(v,digits=30):
    v=Interval.cv(v);assert v.lo>0
    q=10**digits
    l=isqrt(v.lo.numerator*q*q//v.lo.denominator)
    h=isqrt(-((-v.hi.numerator*q*q)//v.hi.denominator))+1
    return Interval(F(l,q),F(h,q))
def log_rat(x,N=80):
    x=F(x);assert x>=1
    z=(x-1)/(x+1);p=z;a=F(0)
    for j in range(N):
        a+=2*p/F(2*j+1);p*=z*z
    return Interval(a,a+2*p/(F(2*N+1)*(1-z*z)))
def outward(x,d=8):
    q=10**d
    a=x.lo.numerator*q//x.lo.denominator
    b=-((-x.hi.numerator*q)//x.hi.denominator)
    def f(t):
        sg='-' if t<0 else '';t=abs(t)
        return sg+str(t//q)+'.'+str(t%q).zfill(d)
    return '['+f(a)+', '+f(b)+']'

sqrt2,sqrt3=sqr(2),sqr(3)
log2,log3=log_rat(2),log_rat(3)
s=4/sqrt3
Ny=Interval(F(22,3))
Np=(44-16*sqrt3)/3
Nm=(44+16*sqrt3)/3
assert F(2)<s.lo<s.hi<F(7,3)

def isprime(n):return n>=2 and all(n%j for j in range(2,isqrt(n)+1))
def prime_base(n):
    for p in range(2,n+1):
        if not isprime(p):continue
        v=n
        while v%p==0:v//=p
        if v==1:return p
    return None
pp={n:prime_base(n) for n in range(2,55) if prime_base(n)}
assert len(pp)==24
base=[(2*j-3,0) for j in range(4)]
def triple(h):return [(2*j-3,h) for j in range(3)]
def reflect(cells):return [(-a,-b) for a,b in cells]
R2,R4=triple(1),triple(2)
L2,L4=reflect(R2),reflect(R4)
centers={'base':base,'R2':R2,'L2':L2,'R4':R4,'L4':L4}
def ratio(x,y):
    assert (x[0]-y[0])%2==0
    return F(3)**((x[0]-y[0])//2)*F(2)**(x[1]-y[1])
exponential_upper=1/(1-ELL)
def hits(A,B):
    out=[]
    for i,x in enumerate(A):
        for j,y in enumerate(B):
            r=ratio(x,y)
            m=max(r,1/r)
            for n in pp:
                difference=max(m/F(n),F(n)/m)
                if difference<exponential_upper:
                    assert difference==1,('near missed prime',A,B,i,j,n,difference)
                    out.append((i,j,n))
    return out
expected={
('base','R2'):[(0,0,2),(1,1,2),(2,2,2)],
('base','R4'):[(0,0,4),(1,1,4),(2,2,4)],
('R2','L2'):[(1,2,4),(2,1,4)],
('R4','L4'):[(1,2,16),(2,1,16)],
('R2','R4'):[(0,0,2),(1,1,2),(2,2,2)],
('R2','L4'):[(1,2,8),(2,1,8)],
('L2','R4'):[(1,2,8),(2,1,8)],
('L2','L4'):[(0,0,2),(1,1,2),(2,2,2)],
}
for (a,b),need in expected.items():assert hits(centers[a],centers[b])==need,(a,b,hits(centers[a],centers[b]))
assert all(ratio(x,y)!=1 for a in centers for b in centers if a!=b for x in centers[a] for y in centers[b])
# Pole null in each triple and the four original cells.
r=S.sqrt(3);s_sym=4/r
X=S.Matrix([[1,0],[-s_sym,1],[1,-s_sym],[0,1]])
L=S.Matrix([[r**j for j in range(4)],[r**(-j) for j in range(4)]])
assert S.simplify(L*X)==S.zeros(2,2)
for h in (1,2):
    moment=S.Matrix([[r**j for j in range(3)], [r**(-j) for j in range(3)]])*S.Matrix([1,-s_sym,1])
    assert S.simplify(moment)==S.zeros(2,1)
# Cross-prime symbols; parity basis gives TWO independent 3x3 matrices.
K2,K4,K8,K16,K3,K9=S.symbols('K2 K4 K8 K16 K3 K9',real=True)
spNy=2+s_sym*s_sym
wp=S.Matrix([1,-s_sym,1]);cp=X*S.Matrix([1,1]);cm=X*S.Matrix([1,-1]);
Np_sym=S.simplify(cp.dot(cp));Nm_sym=S.simplify(cm.dot(cm))
for sign,c,N in [('+',cp,Np_sym),('-',cm,Nm_sym)]:
    D=S.simplify(sum(c[j]*wp[j] for j in range(3)))
    assert S.simplify(D-(2-2*s_sym+s_sym*s_sym if sign=='+' else 2+2*s_sym+s_sym*s_sym))==0
    # q(u_sigma, v_{h,sigma}) = -sqrt(2) K_{2^h} D/sqrt(N Ny).
    # q(v2_sigma,v4_sigma) = -K2 +/- (2s/Ny) K8.
    assert S.simplify(wp.dot(wp)-spNy)==0
assert S.simplify((4*s_sym*K3-2*K9)/spNy-(S.Rational(6,11)*s_sym*K3-S.Rational(3,11)*K9))==0

# Independent 6-mode source Gram and full prime-power matrix check.
# Six modes share 16 individually disjoint, translated unit cells; the
# coefficient map is exact, including all parity and normalization factors.
allcells=base+R2+L2+R4+L4
assert len(allcells)==len(set(allcells))==16
cplus=X*S.Matrix([1,1]);cminus=X*S.Matrix([1,-1])
Ny_sym=spNy
coeffs=[]
for c,N in [(cplus,Np_sym),(cminus,Nm_sym)]:
    vec=S.zeros(16,1)
    for j in range(4):vec[j]=c[j]/S.sqrt(N)
    coeffs.append(vec)
for posR,posL in [(4,7),(10,13)]:
    for sign in (1,-1):
        vec=S.zeros(16,1)
        for j,val in enumerate([1,-s_sym,1]):
            vec[posR+j]=val/S.sqrt(2*Ny_sym)
            vec[posL+j]=sign*val/S.sqrt(2*Ny_sym)
        coeffs.append(vec)
Z=S.Matrix.hstack(*coeffs)
assert S.simplify(Z.T*Z-S.eye(6))==S.zeros(6)
# Gamma diagonal is *not* included here: this is the arithmetic term only.
Qsource=S.zeros(16)
prime_symbols={n:S.Symbol('log'+str(pp[n]),positive=True)/S.sqrt(n) for n in pp}
for i,x in enumerate(allcells):
    for j,y in enumerate(allcells):
        rrat=ratio(x,y); m=max(rrat,1/rrat)
        for n in pp:
            if max(m/F(n),F(n)/m)<exponential_upper:
                assert m==n
                Qsource[i,j]-=prime_symbols[n]
Q6=S.simplify(Z.T*Qsource*Z)
assert S.simplify(Q6-Q6.T)==S.zeros(6)
for i,j in [(0,1),(0,3),(0,5),(1,2),(1,4),(2,3),(2,5),(3,4),(4,5)]:
    assert S.simplify(Q6[i,j])==0,(i,j,Q6[i,j])
for sign,idx,N,D in [('+',0,Np_sym,2-2*s_sym+s_sym**2),('-',1,Nm_sym,2+2*s_sym+s_sym**2)]:
    partner=2 if sign=='+' else 3
    partner4=4 if sign=='+' else 5
    for t,n in [(partner,2),(partner4,4)]:
        lhs=S.simplify(Q6[idx,t]);rhs=S.simplify(-S.sqrt(2)*prime_symbols[n]*D/S.sqrt(N*Ny_sym))
        assert S.simplify(lhs*lhs-rhs*rhs)==0
        assert lhs.is_negative and rhs.is_negative
    assert S.simplify(Q6[partner,partner4]+prime_symbols[2]-(1 if sign=='+' else -1)*2*s_sym*prime_symbols[8]/Ny_sym)==0
    for t,n in [(partner,4),(partner4,16)]:
        assert S.simplify(Q6[t,t]-(S.Rational(7,11)*S.sqrt(3)*prime_symbols[3]+(1 if sign=='+' else -1)*2*s_sym*prime_symbols[n]/Ny_sym))==0

# All translated centers are within physical window I, even for the outermost x2+log4.
# log4 +(1/2)log3 <2 from 4 sqrt(3)<e^2; exact rational series.
exp2_partial=sum((F(2)**j/F(__import__('math').factorial(j)) for j in range(10)),F(0))
assert F(4)*sqrt3.hi<exp2_partial
# For distinct groups, ratio >= 9/8, hence gap t>=log(9/8)-ell:
all_ratios=[max(ratio(x,y),1/ratio(x,y)) for a in centers for b in centers if a!=b for x in centers[a] for y in centers[b]]
assert min(all_ratios)==F(9,8)
# w(t) <= 1/(1-e^{-2t}) < 5, using e^{2ell} <= (1-2ell)^{-1}.
assert 1/(1-F(64,81)/(1-2*ELL))<5
# For base versus 2-cells the ratio>=3/2, so w<2 is safe.
assert min(max(ratio(x,y),1/ratio(x,y)) for x in base for y in R2)==F(3,2)
assert 1/(1-F(4,9)/(1-2*ELL))<2

# Exact prime strengths and conservative Gamma crossing bounds.
kap2=log2/sqrt2;kap4=log2/2;kap8=log2/(2*sqrt2);kap16=log2/4
S1=2+s
assert S1.hi<F(13,3)
probes={}
for tag,N,D,cabs in [('+',Np,2-2*s+s*s,2*s),('-',Nm,2+2*s+s*s,4+2*s)]:
    denom=sqrt_interval(N*Ny)
    for h,kappa,kernel in [(2,kap2,2),(4,kap4,5)]:
        prime=kappa*D/denom
        gamma=kernel*ELL*cabs*S1/denom
        assert prime.lo>gamma.hi
        probes[(tag,h)]=(prime,gamma)
        print(tag,'translated by',h,'prime',outward(prime),'Gamma upper',outward(gamma))

measured={}
for tag in ('+','-'):
    p2,g2=probes[(tag,2)];p4,g4=probes[(tag,4)]
    lo=2*((p2.lo-g2.hi)**2+(p4.lo-g4.hi)**2)
    hi=2*((p2.hi+g2.hi)**2+(p4.hi+g4.hi)**2)
    measured[tag]=Interval(lo,hi)
    print(tag,'certified squared coupling captured by V2+V4:',outward(measured[tag]))
assert measured['+'].lo>F(13,100) and measured['+'].hi<F(7,50)
assert measured['-'].lo>F(29,50) and measured['-'].hi<F(3,5)

# For cross v2-v4, two Gamma correlations (RR and RL), each < 5ell(2+s)^2/Ny.
gamma_24=10*ELL*S1*S1/Ny
prime24_ub=kap2+(2*s/Ny)*kap8
assert (prime24_ub+gamma_24).hi<F(2,3)
# The h=4 diagonal is bounded below exactly as h=2 in T1D-J:
# same Gamma diagonal and internal 3/9, while |q(yR4,yL4)|
# has smaller prime16 than prime4, and Gamma cross <1/1000
# (smallest r4/reflected ratio >=27/16 > 4/3).
assert kap16.hi<kap4.lo
assert (3*ELL*S1*S1/Ny).hi<F(1,1000)
assert F(679,100)+F(7,11)*F(109,100)-F(1,5000)-F(1,4)>F(36,5)

# Conservative *full Weil* six-dimensional lower bound by parity
# and strict Gershgorin dominance: per sector (u_sigma,v2_sigma,v4_sigma),
# diag > (7,36/5,36/5), offdiag < (sqrt2/4 or sqrt2/2,
# sqrt2/5 or 2sqrt2/5, 2/3).
for tag,a,b in [('+',F(1,4),F(1,5)),('-',F(1,2),F(2,5))]:
    assert 7-F(3,2)*(a+b)>F(11,2)
    assert F(36,5)-F(3,2)*a-F(2,3)>F(11,2)
    assert F(36,5)-F(3,2)*b-F(2,3)>F(11,2)
    # stronger exact magnitude checks for u-vh from interval enclosures
    for h,bound in [(2,a),(4,b)]:
        p,g=probes[(tag,h)]
        assert p.hi+g.hi<bound

# Conditional on independent T1D-I two-sided analytic beta^2 results:
# the part of r_sigma outside W=V2+V4 stays quantitatively nonzero.
assert F(16,5)-measured['+'].hi>F(153,50) # >3.06
assert F(12,5)-measured['-'].hi>F(9,5) # >1.8
assert max(measured['+'].hi,measured['-'].hi)/7<F(3,35)
# Using earlier ||R||>320/1393, the full penalty norm cannot be
# attained by compression onto W; this comparison uses T1D-I.
assert F(320,1393)-F(3,35)==F(1003,6965)
print('PASS T1D-K: six-dimensional genuine pole-null physical compression q > 11/2')
print('PASS T1D-K: translated 2/4 prime collisions, mixed 8 and reflected 16')
print('PASS T1D-K: measured squared residual + in (13/100,7/50); - in (29/50,3/5)')
print('CONDITIONAL on T1D-I: squared unmeasured residuals >153/50 and >9/5;')
print('                   ||R||-||Pi_W R Pi_W|| >1003/6965.')
print('OPEN: whole-complement coercivity, conditioned Sonin transport, RH.')
