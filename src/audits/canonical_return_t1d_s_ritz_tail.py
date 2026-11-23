#!/usr/bin/env python3
"""CR-INT-T1D-S: exact certified parity Ritz lower bounds and R.11 tail barrier.

Tests exact original prime-power collision algebra on the authentic T1D-K
smooth complementary modes, then rigorous rational lower quotient bounds
using independently stated T1C/K Gamma estimates. Also verifies the
Dirichlet min-max / conservative Galerkin tail arithmetic.

NO computation of the top compact eigenvalues or proof of whole-complement
positivity, no Lean and no RH.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp

class I:
    def __init__(self,a,b=None):
        self.lo=F(a);self.hi=F(a if b is None else b)
        assert self.lo<=self.hi
    def __add__(self,x):
        x=x if isinstance(x,I) else I(x)
        return I(self.lo+x.lo,self.hi+x.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,x):return self+(-x if isinstance(x,I) else I(-x))
    def __rsub__(self,x):return I(x)-self
    def __mul__(self,x):
        x=x if isinstance(x,I) else I(x)
        q=[a*b for a in (self.lo,self.hi) for b in (x.lo,x.hi)]
        return I(min(q),max(q))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=x if isinstance(x,I) else I(x)
        assert x.lo>0
        return self*I(F(1,x.hi),F(1,x.lo))
    def __rtruediv__(self,x):return I(x)/self
    def __repr__(self):
        d=10**9
        a=(self.lo.numerator*d)//self.lo.denominator
        b=-((-self.hi.numerator*d)//self.hi.denominator)
        return f'[{a/d:.9f},{b/d:.9f}] (rounded display only)'

def root(n, precision=10**30):
    z=isqrt(n*precision*precision)
    assert z*z<=n*precision*precision<(z+1)*(z+1)
    return I(F(z,precision),F(z+1,precision))

def log_rat(q,N=90):
    q=F(q);assert q>=1
    z=(q-1)/(q+1); term=z;sm=F(0)
    for j in range(N):
        sm+=2*term/F(2*j+1);term*=z*z
    return I(sm,sm+2*term/(F(2*N+1)*(1-z*z)))

ln2=log_rat(2,120);ln3=log_rat(3,120)
s=4/root(3);ny=I(F(22,3));a=2*s/ny
k2=ln2/root(2);k4=ln2/2;k8=ln2/root(8);k16=ln2/4

# Exact pairwise prime form from the original 12 disjoint unit cells;
# support centers are encoded x = e3/2 log(3) + e2 log(2).
# For each parity use (v2 + v4)/sqrt(2), with each cell's amplitude
# c_j/(2 sqrt(N_y)); the common factor is accounted for in Q.
z=sp.symbols('z',positive=True)  # log 2
w=sp.symbols('w',positive=True)  # log 3
r=sp.sqrt(3); ss=4/r; NN=sp.Rational(22,3)

def modes(eps):
    cells=[]
    for a2 in [1,2]:
        for reflection in [+1,-1]:
            for j,c in enumerate([1,-ss,1]):
                cells.append(((reflection*(2*j-3),reflection*a2),
                              sp.simplify(eps if reflection<0 else 1)*c))
    assert len(set(x for x,_ in cells))==12
    return cells

def prime_log_weight(p,q):
    dx=p[0]-q[0]
    assert dx%2==0
    dy=p[1]-q[1]
    ratio=F(3)**(dx//2)*F(2)**dy
    if ratio<1:ratio=1/ratio
    if ratio.denominator!=1:return 0
    n=ratio.numerator
    if n<=1:return 0
    for base,logbase in [(2,z),(3,w)]:
        k=n;num=0
        while k%base==0:
            k//=base;num+=1
        if k==1:
            return logbase/sp.sqrt(n)
    return 0

for eps,label,threshold in [(1,'even',F(7,25)),(-1,'odd',F(31,100))]:
    cs=modes(eps)
    raw=sum((sp.simplify(c*d)*prime_log_weight(x,y)
             for i,(x,c) in enumerate(cs) for y,d in cs[i+1:]),sp.Integer(0))
    # Original A quadratic is 2*(sum unordered c_i c_j kappa) / (4N_y).
    A_exact=sp.simplify(raw/(2*NN))
    A_expected= -sp.Rational(7,11)*w - eps*(2*ss/NN)*(z/2+z/4)/2 \
                + z/sp.sqrt(2) - eps*(2*ss/NN)*z/sp.sqrt(8)
    assert sp.simplify(A_exact-A_expected)==0,(label,sp.simplify(A_exact-A_expected))
    Aval=-F(7,11)*ln3 - eps*a*(k4+k16)/2 + k2-eps*a*k8
    numerator=6+Aval  # R>=0, so C>=6+A
    # Source Gamma bound: per diagonal <=g0/J+(4s+2)*ell/Ny+
    # 3 ell (2+s)^2/Ny; the cross gamma between v2,v4 is
    # <=10 ell (2+s)^2/Ny.  For their unit average:
    ell=F(1,10000)
    gamma_extra=ell*((4*s+2)+13*(2+s)*(2+s))/ny
    assert gamma_extra.hi<F(1,100)
    # T1C explicitly gives g0/J < 13.44129808 < 13.442;
    # hence d=6+gamma[w] < 6+13.442+0.01 <39/2.
    assert F(6)+F(13442,1000)+F(1,100)<F(39,2)
    assert numerator.lo>threshold*F(39,2),(label,numerator.lo,threshold*F(39,2))
    print(f'PASS T1D-S {label}: exact 12-cell original prime identity; 6+<Aw,w> = {numerator}')
    print(f'PASS T1D-S {label}: Ritz Rayleigh >{threshold} using 0<=R and d<39/2')

# Galerkin estimator in T1D-R: E_N=30 sqrt(8/(5d_{N+1}))+15/d_{N+1}.
# If this specific upper estimator is <=1/2, its FIRST positive term
# alone requires d_{N+1} >= 5760. From the first N+5 Dirichlet modes,
# by codimension <=4 and minmax, d_{N+1} <=46+log(1+pi(N+5)/4).
assert F(30)**2 *F(8,5) /F(1,2)**2==5760
assert 5760-46==5714
# No literal huge transcendental output is treated as a computable N.
# The necessary rank inequality is log(1+pi(N+5)/4) >=5714.
for N in [10**6,10**8,10**12]:
    # pi<22/7 and 1+(pi/4)*(N+5) < 3N for N >=10.
    assert N>=10 and 1+F(22,7)*F(N+5,4) <3*N
    # log(3N) <= log 3 + k log 10; write N=10^q here.
    q=len(str(N))-1
    assert N==10**q
    log10=3*ln2+log_rat(F(5,4),100)
    log3N=ln3+q*log10
    assert 46+log3N.hi<5760
    print(f'PASS T1D-S: N=10^{q} cannot make the stated R.11 tail estimator <=1/2')

print('NECESSARY for R.11 estimator <=1/2: N >= (4/pi)*(exp(5714)-1)-5.')
print('OPEN: actual leading compact spectrum, sharper operator-adapted tail, full Schur positivity, RH.')
