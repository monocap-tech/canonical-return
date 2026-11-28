#!/usr/bin/env python3
"""CR-ALP-029: rigorous Fraction-interval POSITIVE SIX-MODE lower certificate.

The original signed Weil form on the physical pole-null w_3,...,w_8 is
Gamma g(0)*M + sum_{k>=0} J_{beta_k} + unchanged seven paired primes.
Each J_beta is PSD. We strictly verify H_K := g(0)M+sum_{k<K}J_beta+P
is positive definite via an OUTWARD-RATIONAL LDL^T. Therefore the complete
six-mode original signed form is positive definite by the analytic identity.

Exact checks employ the source-attached CR-ALP-028 Fraction interval code.
Classical analytic inputs: digamma exponential kernel, sine Laplace integral,
Euler-Maclaurin alternating harmonic bound and exact original EF. The Python
script audits rational enclosures, not those analytical theorems or Lean.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import importlib.util
import json
import sys

HERE=Path(__file__).resolve().parent
PRIME_AUDIT=HERE.parent/'cr_alp_028'/'cr_alp_028_prime_interval_and_gram_audit.py'
if not PRIME_AUDIT.exists():
    PRIME_AUDIT=HERE/'cr_alp_028_prime_interval_and_gram_audit.py'
spec=importlib.util.spec_from_file_location('orig_prime', PRIME_AUDIT)
prime=importlib.util.module_from_spec(spec)
spec.loader.exec_module(prime)
I=prime.I
A=prime.A
PI=prime.pi

def log_interval_pos(z, terms=65):
    """Certified 2*atanh((z-1)/(z+1)), with positive 1<=z<=2."""
    z=I(z)
    assert 1<=z.lo<=z.hi<=2
    u=(z-I(1))/(z+I(1))
    # For positive u, lower partial and geometrically bounded upper remainder.
    lo=I(0); powu=u
    for j in range(terms):
        lo=lo+2*powu/I(2*j+1)
        powu=powu*u*u
    tail=2*powu/(I(2*terms+1)*(I(1)-u*u))
    return I(lo.lo,(lo+tail).hi)

def exp_rational_pos(z=F(11,10),terms=75):
    assert z>0 and F(z)<2
    s=F(0)
    for k in range(terms+1):
        s+=z**k/F(factorial(k))
    first=z**(terms+1)/factorial(terms+1)
    tail=first/(1-z/F(terms+2))
    return I(s,s+tail)

# Euler--Maclaurin, classical strict alternating bounds, n>=1:
# H_n - ln(n) - [1/(2n)-1/(12n²)+1/(120n⁴)] < gamma
# < H_n - ln(n) - [1/(2n)-1/(12n²)+1/(120n⁴)-1/(252n⁶)]
# (Here the left and right ordering derives from D=H-ln n-gamma
#  being between D6 and D4.)
n=256
Hn=sum((F(1,k) for k in range(1,n+1)), F(0))
lnn=prime.log_int(n)
D4=F(1,2*n)-F(1,12*n**2)+F(1,120*n**4)
D6=D4-F(1,252*n**6)
GAMMA_EULER=I((I(Hn)-lnn-I(D4)).lo,(I(Hn)-lnn-I(D6)).hi)
assert GAMMA_EULER.lo>F(57,100) and GAMMA_EULER.hi<F(58,100)
ln_pi=prime.log_int(2)+log_interval_pos(PI/I(2))
q=I(1)/exp_rational_pos()
assert q.lo>F(3,10) and q.hi<F(2,5)
g0=-GAMMA_EULER-PI/I(2)-I(3)*prime.log_int(2)-ln_pi
assert g0.lo>F(-6) and g0.hi<F(-5)
print('PASS: Euler constant & log(pi) enclosed by exact rational series')
print('g(0) in',str(float(g0.lo)),str(float(g0.hi)))

# Convert all source constants to outward-rounded fixed-point INTEGER intervals.
# Each integer denotes an exact multiple of 10^-42; all operations are
# exact integer additions and rigorous floor/ceiling divisions.  This
# avoids hundreds of thousands of costly Fraction gcd normalisations.
SCALE=10**42
class FI:
    __slots__=('lo','hi')
    @staticmethod
    def make(lo,hi):
        v=object.__new__(FI);v.lo,v.hi=lo,hi;return v
    def __init__(self,x=0):
        if isinstance(x,FI):self.lo,self.hi=x.lo,x.hi;return
        if isinstance(x,I):
            self.lo=(x.lo.numerator*SCALE)//x.lo.denominator
            self.hi=-((-x.hi.numerator*SCALE)//x.hi.denominator)
            return
        x=F(x)
        self.lo=(x.numerator*SCALE)//x.denominator
        self.hi=-((-x.numerator*SCALE)//x.denominator)
    def __add__(self,b):
        b=FI(b);return FI.make(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return FI.make(-self.hi,-self.lo)
    def __sub__(self,b):return self+-FI(b)
    def __rsub__(self,b):return FI(b)+-self
    def __mul__(self,b):
        b=FI(b);p=(self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi)
        a,c=min(p),max(p)
        return FI.make(a//SCALE,-((-c)//SCALE))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=FI(b)
        assert b.lo>0, 'nonpositive denominator interval'
        lows=[(a*SCALE)//d for a in (self.lo,self.hi) for d in (b.lo,b.hi)]
        highs=[-((-a*SCALE)//d) for a in (self.lo,self.hi) for d in (b.lo,b.hi)]
        return FI.make(min(lows),max(highs))
    def __repr__(self):return f'[ {self.lo}/{SCALE}, {self.hi}/{SCALE} ]'

JS=(3,4,5,6,7,8)
LAM={j:FI(prime.lam(j)) for j in JS}
NU={j:FI(prime.nu[j]) for j in (1,2,*JS)}
q=FI(q)
g0=FI(g0)
A=FI(A)
M={}
for i in JS:
 for j in JS:
    M[(i,j)]=FI(int(i==j))+(LAM[i]*LAM[j] if i%2==j%2 else FI(0))

H={}
for i in JS:
 for j in JS:
    if i%2!=j%2: H[(i,j)]=FI(0)
    else: H[(i,j)]=g0*M[(i,j)]+FI(prime.Pw(i,j))

K=int(sys.argv[1]) if len(sys.argv)>1 else 30000
assert 1<=K<=300000
q4=q*q*q*q
qpower=q
for k in range(K):
    beta=F(4*k+1,2)
    be=FI(beta)
    D={j:FI(beta*beta)+NU[j]*NU[j] for j in NU}
    diag={j:FI(2)*NU[j]*NU[j]/(be*D[j]) for j in NU}
    h={j:NU[j]/D[j] for j in NU}
    for js,anchor,eps in [((3,5,7),1,-1),((4,6,8),2,1)]:
        factor=FI(2)/A*(FI(1)-FI(eps)*qpower)
        t={j:h[j]-LAM[j]*h[anchor] for j in js}
        for ii,i in enumerate(js):
            for j in js[:ii+1]:
                value=(diag[i] if i==j else FI(0))+LAM[i]*LAM[j]*diag[anchor]-factor*t[i]*t[j]
                H[(i,j)]=H[(i,j)]+value
                if i!=j:H[(j,i)]=H[(i,j)]
    qpower=qpower*q4
    if (k+1)%10000==0:print('summed',k+1,'strictly positive translation-difference kernels',flush=True)

# Outward-rational LDL^T in the exact two parity blocks. All pivots
# must have strictly positive rational lower endpoints.
def ldl_pivots(js):
    L={}
    piv=[]
    for i,ji in enumerate(js):
        di=H[(ji,ji)]
        for h in range(i):
            di=di-L[(i,h)]*L[(i,h)]*piv[h]
        assert di.lo>0,(K,js,i,str(di))
        piv.append(di)
        for j in range(i+1,len(js)):
            numerator=H[(js[j],ji)]
            for h in range(i):
                numerator=numerator-L[(j,h)]*L[(i,h)]*piv[h]
            L[(j,i)]=numerator/di
    return piv,L

pivot_odd,Lodd=ldl_pivots((3,5,7))
pivot_even,Leven=ldl_pivots((4,6,8))

def gen_lower_bound(js,pivs,L):
    # Unit-lower-triangular LDL: ||L^{-1}||_F^2 bounds L^{-T}.
    # For 3x3, the inverse has (2,0)=l21*l10-l20.
    ell10,ell20,ell21=L[(1,0)],L[(2,0)],L[(2,1)]
    frob=FI(3)+ell10*ell10+ell21*ell21+(ell21*ell10-ell20)*(ell21*ell10-ell20)
    trace_mass=sum((M[(j,j)] for j in js),FI(0))
    dmin=min(p.lo for p in pivs)
    # Since M PSD, M <= tr(M) I. Hence H>=dmin/frob I
    # and H>=lambda M for lambda=dmin/(frob*trM).
    return F(dmin,SCALE)/(F(frob.hi,SCALE)*F(trace_mass.hi,SCALE))

lower_odd=gen_lower_bound((3,5,7),pivot_odd,Lodd)
lower_even=gen_lower_bound((4,6,8),pivot_even,Leven)
lower_generalized=min(lower_odd,lower_even)
print("CERTIFIED generalized eigenvalue lower bound: >",format(float(lower_generalized),".12g"))
assert lower_generalized>F(5,10**9),"Unable to certify 5e-9 generalized margin"
print('CERTIFIED strictly positive 6x6 original signed partial-Gamma + full-primes matrix')
for parity,pivs in [('odd-index/even-physical',pivot_odd),('even-index/odd-physical',pivot_even)]:
 print(parity,':')
 for p in pivs:print('  pivot lower:',format(p.lo/SCALE,'.14g'),'upper:',format(p.hi/SCALE,'.14g'))
assert all(p.lo>0 for p in pivot_odd+pivot_even)

# Original six-mode H_full = H_K + sum_{k>=K} J_beta, with EACH J_beta PSD,
# so certified H_K >> 0 gives H_full >> 0. No tail or global W claim.
report={
 'stage':'CR-ALP-029',
 'matrix':'original signed full Weil on w_3 through w_8, strictly positive',
 'source':'actual seven paired von Mangoldt translations + complete digamma gamma',
 'finite_gamma_terms':K,
 'precision_denominator':SCALE,
 'certified_generalized_eigenvalue_lower_bound_fraction':str(lower_generalized),
 'certified_generalized_eigenvalue_lower_bound_gt_5e_9':True,
 'positive_LDL_lower_pivots':{
  'index_odd':[str(F(p.lo,SCALE)) for p in pivot_odd],
  'index_even':[str(F(p.lo,SCALE)) for p in pivot_even],
 },
 'analytic_input':'digamma positive kernel expansion, sine Laplace formula, classical Euler--Maclaurin gamma constant',
 'status':'EXACT OUTWARD-FRACTION ARITHMETIC PASS + ANALYTIC INPUT; NOT LEAN; NOT FULL COLLAR POSITIVITY',
}
(HERE/'CR-ALP-029_certified_six_mode_result.json').write_text(json.dumps(report,indent=2)+'\n')
print('CR-ALP-029 6x6 ORIGINAL SIGNED WEIL POSITIVITY INTERVAL AUDIT PASS')
