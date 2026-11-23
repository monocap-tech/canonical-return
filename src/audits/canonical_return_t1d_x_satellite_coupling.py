#!/usr/bin/env python3
"""CR-INT-T1D-X: exact rational certificate for central/satellite mixed Weil
pairings, original complete prime-power census and positive two-mode Gram.

The fixed T1D-W mollifier is arbitrary smooth even nonnegative probability
supported in [-1e-7, 1e-7]; the inequalities hold uniformly for that class.
Gamma diagonal brackets are explicitly IMPORTED from T1D-W and should be
rechecked with its independent program. Gamma OFF-diagonal bounds are new
analytical integration-by-parts calculations with exact rational checks here.
The additional q-R coupling bound IMPORTS T1D-I's residual estimate.
Not Lean, no full-complement sign, RH or completed Suzuki-Connes transport.
"""
from fractions import Fraction as F
from math import isqrt, factorial
import sympy as sp

class Interval:
    def __init__(self,lo,hi=None):
        self.lo=F(lo);self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,other):
        other=other if isinstance(other,Interval) else Interval(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,other):return self+(-other if isinstance(other,Interval) else Interval(-other))
    def __rsub__(self,other):return Interval(other)-self
    def __mul__(self,other):
        other=other if isinstance(other,Interval) else Interval(other)
        choices=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(choices),max(choices))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=other if isinstance(other,Interval) else Interval(other)
        assert other.lo>0
        return self*Interval(1/other.hi,1/other.lo)
    def __repr__(self):
        return f'({float(self.lo):.12f},{float(self.hi):.12f}) diagnostic'

def log_rat(r,N=90):
    r=F(r); assert r>=1
    h=(r-1)/(r+1);hh=h*h;term=h;tot=F(0)
    for j in range(N):
        tot+=2*term/F(2*j+1);term*=hh
    return Interval(tot,tot+2*term/(F(2*N+1)*(1-hh)))

def root(n,den=10**24):
    lo=isqrt(n*den*den)
    assert lo*lo<=n*den*den<(lo+1)*(lo+1)
    return Interval(F(lo,den),F(lo+1,den))

def prime_base(n):
    for p in range(2,n+1):
        if n%p:continue
        if any(p%q==0 for q in range(2,isqrt(p)+1)):continue
        m=n
        while m%p==0:m//=p
        if m==1:return p
    return None

pp={n:prime_base(n) for n in range(2,55) if prime_base(n)}
assert len(pp)==24 and list(pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]
kappa={n:log_rat(p)/root(n) for n,p in pp.items()}
# All actual cross-shifts must remain inside the original I=(-2,2).
# The central interval [-a-eps, a+eps], a=1/4, and satellite copies
# at +/- log3, are strictly inside I and disjoint from E (at +/-(log3)/2,
# +/-(3log3)/2 with narrow carrier radius <=1/20000).
eps=F(1,10**7);a=F(1,4);L=2*(a+eps)
ln3=log_rat(3)
assert ln3.lo>F(1098,1000)
assert ln3.lo+a*(-1)-eps>F(4,5) # satellite left edge >0.8
assert ln3.hi+a+eps<F(3,2)      # satellite right edge <1.5
assert ln3.lo-a-eps>log_rat(3).hi/2+F(1,5) # away from selected inner cell
assert 3*ln3.lo/2-F(1,20000)>ln3.hi+a+eps  # away from selected outer cell
assert L<log_rat(2).lo
# Actual central-to-positive-satellite prime collisions:
# |log n - log3| < L exactly when n=2,3,4.
central=[]
for n in pp:
    if n==3: central.append(n); continue
    h=log_rat(F(max(n,3),min(n,3)))
    if h.hi<L:central.append(n)
    else:assert h.lo>L
assert central==[2,3,4]
# Satellite-wing pair separation is 2 log3=log9, so original prime
# overlap when |log n-log9|<L exactly when n=7,8,9,11,13.
wing=[]
for n in pp:
    if n==9:wing.append(n);continue
    h=log_rat(F(max(n,9),min(n,9)))
    if h.hi<L:wing.append(n)
    else:assert h.lo>L
assert wing==[7,8,9,11,13]

# Autocorrelation of the UNMOLLIFIED compact polynomial f=(d²-1/4)b.
# For shift 0<=h<1/2 the exact overlap is a polynomial in h:
# R_f(h)=integral_{h-a}^a f(x) f(x-h) dx. The symmetrical change of
# variable t=x-h/2 reduces it to even powers integrated exactly.
x,t,h=sp.symbols('x t h', real=True)

def toF(q):
    q=sp.cancel(q); assert q.is_Rational, q
    return F(int(q.p),int(q.q))

def raw_data(parity):
    b=sp.expand(x**parity*(1-16*x*x)**3)
    f=sp.expand(sp.diff(b,x,2)-b/4)
    P=sp.Poly(sp.expand(f.subs(x,t+h/2)*f.subs(x,t-h/2)),t)
    Rh=sp.Poly(sp.expand(sum(c*2*(sp.Rational(1,4)-h/2)**(i+1)/(i+1)
              for (i,),c in P.terms() if i%2==0)),h)
    m=toF(Rh.eval(0))
    grad=toF(sp.integrate(sp.diff(f,x)**2,(x,-sp.Rational(1,4),sp.Rational(1,4))))
    l1b=toF(sp.integrate((1-16*x*x)**3 if parity==0 else x*(1-16*x*x)**3,
                      (x,0,sp.Rational(1,4)))) * (2 if parity==0 else 2)
    # parity=0 is even, parity=1 odd so 2* integral_0^a |b| exact
    return Rh,m,grad,l1b

def eval_poly_interval(poly,q):
    # Exact interval-Horner: rational polynomial on a rigorously
    # enclosed rational logarithmic shift.
    v=Interval(0)
    for c in poly.all_coeffs():
        v=v*q+Interval(toF(c))
    return v

# Bounds on the fourth-order Gamma kernel H(s)=(d²-1/4)^2 w(s).
# H(s)=e^{-s/2} * [16 S4(q)+16 S3(q)+4 S2(q)], q=e^{-2s}.
# For near separation s=log3-L > 0.59, q<31/100,
# for far s=2log3-L >1.69, q<1/25.
def exp_lower(z, N=24):
    return sum((F(z)**k/F(factorial(k)) for k in range(N+1)),F(0))
assert ln3.lo-L>F(59,100)
assert 2*ln3.lo-L>F(169,100)
assert exp_lower(F(118,100))>F(100,31)
assert exp_lower(F(338,100))>25

def H_series_majorant(q):
    q=F(q)
    S2=q*(1+q)/(1-q)**3
    S3=q*(1+4*q+q*q)/(1-q)**4
    S4=q*(1+11*q+11*q*q+q**3)/(1-q)**5
    return 16*S4+16*S3+4*S2
assert H_series_majorant(F(31,100))<231
assert H_series_majorant(F(1,25))<3

for parity in (0,1):
    label='even' if parity==0 else 'odd'
    R,m,k,l1b=raw_data(parity)
    if parity==0:
        assert m==F(28194976,15015) and k==F(259883264,385) and l1b==F(8,35)
        derivative_max=19
    else:
        assert m==F(184586,3465) and k==F(23448144,715) and l1b==F(1,64)
        derivative_max=25
    assert k<derivative_max**2*m
    delta=2*eps*derivative_max+eps**2*derivative_max**2
    assert delta<F(1,100000) # physical mollification mass correction
    mass=Interval(1-delta,1+delta)
    def autocorr(r):
        if r==1:return Interval(1)
        shift=log_rat(r)
        assert shift.hi<F(1,2)
        raw=eval_poly_interval(R,shift)/Interval(m)
        return (raw+Interval(-delta,delta))/mass
    # All original cross entries, normalized by the PHYSICAL SMOOTH MASS:
    # p = -q_prime(y,Tlog3 y)/mass, gamma may add independently.
    pc=kappa[3]+kappa[2]*autocorr(F(3,2))+kappa[4]*autocorr(F(4,3))
    # wing beta = q_prime(Tlog3 y,T-log3 y)/mass,
    # prime n=9 with coefficient -kappa9, and n7,8,11,13.
    pw=Interval(0)
    for n in wing:
        corr=Interval(1) if n==9 else autocorr(F(max(n,9),min(n,9)))
        pw-=kappa[n]*corr
    if parity==0:
        assert F(67,100)<pc.lo<pc.hi<F(68,100)
        assert F(27,100)<pw.lo<pw.hi<F(29,100)
        g_lo,g_hi=F(89,100),F(91,100) # independent T1D-W Gamma diagonal intervals
        Rsmall=F(4,1575) # true Schur penalty on central Y, from T1D-W
    else:
        assert F(70,100)<pc.lo<pc.hi<F(71,100)
        assert -F(11,100)<pw.lo<pw.hi<-F(9,100)
        g_lo,g_hi=F(118,100),F(120,100)
        Rsmall=F(9,4375)
    # Exact Gamma OFF bound from integrating L=d²-1/4 twice against w.
    # For B=b*rho, ||B||1<=||b||1 and ||y||²>=m(1-delta).
    near=F(231)*l1b*l1b/(m*(1-delta))
    far=F(3)*l1b*l1b/(m*(1-delta))
    assert near<F(1,100) and far<F(1,1000)
    # Strict conservative complete signed Weil 2x2 Gram in physical
    # orthonormal normalized basis (y,z), z=(T_d y +T_-d y)/sqrt2.
    # q(y,y)=g; q(z,z)=g+gamma_far+pw;
    # q(y,z)=sqrt(2)*(gamma_near-pc).
    # Use sqrt(2) ONLY through exact squares of conservative bounds.
    p_lo,p_hi=(F(67,100),F(68,100)) if parity==0 else (F(70,100),F(71,100))
    w_lo=F(27,100) if parity==0 else -F(11,100)
    diag_z_lower=g_lo+w_lo-F(1,1000)
    cross_lower_factor=p_lo-F(1,100)
    cross_upper_factor=p_hi+F(1,100)
    det_lower=g_lo*diag_z_lower-2*cross_upper_factor**2
    if parity==0:
        assert det_lower>F(7,100)
    else:
        assert det_lower>F(1,5)
    # Combined ORIGINAL source finite Gram strictly positive.
    # The original Schur mixed coupling remains large: use R Cauchy and
    # optional T1D-I global ||R||<155/182. Since sqrt(2)>7/5,
    # sqrt(Rsmall*155/182)<1/20, the modulus is large. True values
    # are even stronger; reported fractions are conservative.
    assert 2>F(7,5)**2
    assert Rsmall*F(155,182)<F(1,20)**2
    schur_cross_lower=F(7,5)*cross_lower_factor-F(1,20)
    assert schur_cross_lower> (F(17,20) if parity==0 else F(9,10))
    # Necessary (not sufficient) satellite Schur diagonal if the entire
    # complement were nonnegative: s[z] >= |s[y,z]|^2 / s[y].
    # Since s[y] < g_hi and |s[y,z]| > lower, this requires
    # s[z] > lower^2 / g_hi. This is CONDITIONAL on T1D-I.
    necessary=schur_cross_lower**2/g_hi
    assert necessary > (F(83,100) if parity==0 else F(69,100))
    assert necessary == (F(190969,227500) if parity==0 else F(52441,75000))
    print('NECESSARY if full Schur positive',label,'satellite s[z] >',necessary,'(conditional on I)')
    print('PASS',label,'raw mass',m,'derivative',k,'||b||1',l1b)
    print('CERTIFIED',label,'original prime center->satellite magnitude p',pc)
    print('CERTIFIED',label,'original prime wing<->wing signed self beta',pw)
    print('CERTIFIED',label,'Gamma near abs/mass <',float(near),'Gamma far abs/mass <',float(far))
    print('CERTIFIED',label,'2x2 Weil Gram determinant >',det_lower,'(positive)')
    print('CONDITIONAL on T1D-I',label,'|(q-R)(y,z)|/(||y||||z||) >',schur_cross_lower)

print('PASS: prime census central [2,3,4], wing [7,8,9,11,13], original 24-prime-power source')
print('SOURCE IMPORT: T1D-W gamma diagonal intervals and central selected Schur bounds; rerun W checker')
print('CONDITIONAL IMPORT: T1D-I global ||R||<155/182 for q-R cross LOWER bounds only')
print('OPEN: full W remainder sign; projector eta_N; whole F Schur positivity; Suzuki-Connes transfer; RH')
