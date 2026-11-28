#!/usr/bin/env python3
"""CR-INT-T1D-W: rational interval Gamma Grams for a fixed tiny-mollifier
central source, exact prime-2/5 selected-central census, and Schur norm bound.

Mathematical inputs stated in the note: Suzuki's Gamma kernel; digamma/Jensen
energy estimate from T1D-U; Schwarz and smooth-mollifier approximation.
All finite interval operations here use integers/Fraction and directed rounding.
No numerical FFT, eigenvalue-upper certificate, full complement theorem or RH.
"""
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import sympy as sp

S=10**31

def fl(q):
    q=F(q)
    return (q.numerator*S)//q.denominator

def ce(q):
    return -fl(-F(q))

class Iv:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None,raw=False):
        if raw:
            self.lo=int(lo);self.hi=int(hi)
        else:
            self.lo=fl(F(lo));self.hi=ce(F(lo if hi is None else hi))
        assert self.lo<=self.hi
    def __add__(self,o):
        o=as_iv(o);return Iv(self.lo+o.lo,self.hi+o.hi,True)
    __radd__=__add__
    def __neg__(self):return Iv(-self.hi,-self.lo,True)
    def __sub__(self,o):return self+-as_iv(o)
    def __rsub__(self,o):return as_iv(o)-self
    def __mul__(self,o):
        o=as_iv(o);z=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return Iv(min(z)//S,-((-max(z))//S),True)
    __rmul__=__mul__
    def scale(self,q):
        q=F(q)
        if q>=0:return Iv(fl(F(self.lo,S)*q),ce(F(self.hi,S)*q),True)
        return Iv(fl(F(self.hi,S)*q),ce(F(self.lo,S)*q),True)
    def __truediv__(self,o):
        if isinstance(o,(Iv,)):
            assert o.lo>0
            return self*Iv(F(S,o.hi),F(S,o.lo))
        return self.scale(1/F(o))
    def __rtruediv__(self,o):return as_iv(o)/self
    def bounds(self):return F(self.lo,S),F(self.hi,S)
    def fmt(self,n=7):
        q=10**n
        lo=(self.lo*q)//S;hi=-((-self.hi*q)//S)
        def f(v):return ('-' if v<0 else '')+str(abs(v)//q)+'.'+str(abs(v)%q).zfill(n)
        return '['+f(lo)+', '+f(hi)+']'

def as_iv(o):return o if isinstance(o,Iv) else Iv(o)

def atan_q(q,N=36):
    q=F(q);assert 0<q<1
    s=sum(((-1)**i*q**(2*i+1)/F(2*i+1) for i in range(N)),F(0))
    t=q**(2*N+1)/F(2*N+1)
    return Iv(s,s+t) if N%2==0 else Iv(s-t,s)

def log_q(q,N=110):
    q=F(q);assert 1<=q<=4**2
    t=(q-1)/(q+1);z=t;tot=F(0)
    for j in range(N):
        tot+=2*z/F(2*j+1);z*=t*t
    return Iv(tot,tot+2*z/(F(2*N+1)*(1-t*t)))

LOG2=log_q(2)
def log_any(q):
    q=F(q);assert q>=1
    m=0
    while q>=2:
        q/=2;m+=1
    return m*LOG2+log_q(q)

def exp_neg_small(x):
    x=F(x);assert 0<=x<=1
    s=F(0);term=F(1)
    for j in range(102):
        s+=term
        term=term*(-x)/F(j+1)
    # 102 terms (j=0..101): last negative -> lower, next positive -> upper
    return Iv(s,s+term)

E1=exp_neg_small(1);EHALF=exp_neg_small(F(1,2));EQUARTER=exp_neg_small(F(1,4))
assert E1.lo>0 and EHALF.lo>0
pi=16*atan_q(F(1,5))-4*atan_q(F(1,239),24)
assert F(31415,10000)<pi.bounds()[0] and pi.bounds()[1]<F(31417,10000)
log4pi=Iv(log_q((4*pi).bounds()[0]).bounds()[0],log_q((4*pi).bounds()[1]).bounds()[1])
H=sum((F(1,j) for j in range(1,1025)),F(0))
euler=Iv(H)-10*LOG2-Iv(F(1,2*1025),F(1,2*1024))
C_GAMMA=log4pi+euler
ratio=(Iv(1)+EHALF)/(Iv(1)-EHALF)
logcoth=Iv(log_q(ratio.bounds()[0]).bounds()[0],log_q(ratio.bounds()[1]).bounds()[1])
BASE=logcoth-C_GAMMA

# Source polynomials and their exact real physical mass, derivative, and correlations.
x,t=sp.symbols('x t',real=True)
a=sp.Rational(1,4)
entries=[]
for parity in [0,1]:
    profiles=[]
    for k in (3,4):
        b=x**parity*(1-16*x*x)**k
        profiles.append(sp.expand(sp.diff(b,x,2)-b/4))
    M=[[sp.integrate(profiles[i]*profiles[j],(x,-a,a)) for j in range(2)] for i in range(2)]
    K=[[sp.integrate(sp.diff(profiles[i],x)*sp.diff(profiles[j],x),(x,-a,a)) for j in range(2)] for i in range(2)]
    assert sp.Matrix(M).det()>0
    for i,j in [(0,0),(0,1),(1,1)]:
        P=sp.Poly(sp.expand(profiles[i]*profiles[j].subs(x,x+t)),x)
        corr=sp.Poly(sp.expand(sum((c/F(int(power[0]+1)))*((a-t)**(power[0]+1)-(-a)**(power[0]+1)) for power,c in P.terms())),t)
        assert corr.eval(0)==M[i][j] and corr.coeff_monomial(t)==0
        coeff=[F(int(sp.numer(corr.nth(k))),int(sp.denom(corr.nth(k)))) for k in range(corr.degree()+1)]
        entries.append((parity,i,j,F(M[i][j]),F(K[i][i]),F(K[j][j]),coeff))

# The Gamma scalar kernel on the physical interval is represented by
# Gamma_ij = (log coth(L/2)-C_Gamma)*M_ij +
#  2 * sum_{m>=0} int_0^L [M_ij exp(-(2m+1)t)
#                            -R_ij(t) exp(-(2m+1/2)t)]dt, L=1/2.
# The tail uses |M_ij-R_ij(t)|<= t^2/2 ||f_i'|| ||f_j'||,
# holding for the equal-parity real zero-extended profiles.
# Hence |tail| <= M/(2alpha)+(2M+D)/(2alpha^2)+2D/alpha^3,
# alpha=2N+1/2, D<= (Kii+Kjj)/2.
N=2048
L=F(1,2)
EPOWER=Iv(1)
GAMMA_RAW={}
for m in range(N):
    # Endpoints at L=1/2: (2m+1)L=m+1/2; (2m+1/2)L=m+1/4.
    Ebeta=EPOWER*EHALF
    Ealpha=EPOWER*EQUARTER
    beta=F(2*m+1)
    alpha=F(4*m+1,2)
    J0=(Iv(1)-Ebeta)/beta
    I=[(Iv(1)-Ealpha)/alpha]
    for k in range(1,20):
        I.append((I[-1]*k-Ealpha.scale(L**k))/alpha)
    for parity,i,j,Mass,Kii,Kjj,cs in entries:
        key=(parity,i,j)
        term=J0*Mass
        for k,c in enumerate(cs):
            if c:term=term-I[k]*c
        GAMMA_RAW[key]=GAMMA_RAW.get(key,Iv(0))+2*term
    EPOWER=EPOWER*E1

for parity,i,j,Mass,Kii,Kjj,cs in entries:
    key=(parity,i,j)
    D=(Kii+Kjj)/2
    alpha=F(4*N+1,2)
    tail=Mass/(2*alpha)+(2*Mass+D)/(2*alpha*alpha)+2*D/(alpha**3)
    GAMMA_RAW[key]=BASE*Mass+GAMMA_RAW[key]+Iv(-tail,tail)

# For a fully specified mollifier rho*(x)=Z^-1 exp(-1/(1-(x/e)^2))
# for |x|<e with e=10^-7, any shape-independent even unit probability
# with this support satisfies the following certified Gamma/mass bounds.
eps=F(1,10**7)

def gamma_jensen_upper(V):
    V=F(V);q=F(1,4)
    terms=sum((V/((n+q)*(4*(n+q)**2+V)) for n in range(51)),F(0))
    return Iv(-F(5372,1000)+terms)+log_any(1+V/(4*(50+q)**2)).scale(F(1,2))

# Source lemma: d[f]=gamma[f]+6||f||^2, and m_Gamma+6>=0.
# For delta_i=(rho*f_i)-f_i, ||delta_i||<=eps||f_i'||,
# ||delta_i'||<=2||f_i'||. The jointly increasing perspective
# M*(6+gamma_Jensen(K/M)) gives d[delta_i] <= 30 eps^2 Kii.
assert (gamma_jensen_upper(F(4)/eps**2)+6).hi<30*S
assert (gamma_jensen_upper(F(45)**2)+6).hi<12*S

def sqrt_up(q):
    q=F(q);assert q>=0
    z=isqrt(q.numerator*S*S//q.denominator)
    return F(z+1,S)

def sqrt_lo(q):
    q=F(q);assert q>=0
    z=isqrt(q.numerator*S*S//q.denominator)
    return F(z,S)

G_MOL={}
M_MOL={}
print('PASS source: original polynomial mass Gram and exact real correlation polynomials (6 entries)')
print('PASS source: exact rational Gamma kernel series 2048 terms; analytic tail bounded with K derivatives')
for parity,i,j,Mass,Kii,Kjj,cs in entries:
    key=(parity,i,j)
    mi=entries[[e[:3] for e in entries].index((parity,i,i))][3]
    mj=entries[[e[:3] for e in entries].index((parity,j,j))][3]
    norml=sqrt_lo(mi*mj); normu=sqrt_up(mi*mj)
    assert norml>0
    gm=GAMMA_RAW[key]
    # Gamma-mollification bound through shifted Gamma-form Cauchy-Schwarz:
    # d_delta <= 30 eps² K, d_raw <=12 M,
    # |d(rho*f_i,rho*f_j)-d(f_i,f_j)|
    # <= sqrt(d_delta_i*d_raw_j)+sqrt(d_delta_j*d_raw_i)+sqrt(d_delta_i*d_delta_j).
    v1=sqrt_up(30*eps*eps*Kii*12*mj)
    v2=sqrt_up(30*eps*eps*Kjj*12*mi)
    v3=sqrt_up(30*eps*eps*Kii*30*eps*eps*Kjj)
    # Physical mass difference via ||delta_i|| <= eps sqrt(Kii).
    dm=eps*(sqrt_up(Kii*mj)+sqrt_up(Kjj*mi))+eps*eps*sqrt_up(Kii*Kjj)
    deviation=v1+v2+v3+6*dm
    gmol=gm+Iv(-deviation,deviation)
    G_MOL[key]=gmol
    M_MOL[key]=Iv(Mass-dm,Mass+dm)
    # normalized by RAW masses (fixed exact); physical smoothed masses differ
    # at relative scale <1e-5, stated separately in note.
    lo,hi=gmol.bounds();q=Iv(lo/norml,hi/norml) if lo>=0 else Iv(lo/norml,hi/norml)
    # correct division by positive interval: evaluate all four corners
    divnorm=Iv(norml,normu)
    q=gmol/divnorm
    print(f'CERTIFIED gamma_{{{("even" if parity==0 else "odd")},{i+3}{j+3}}}/sqrt(M0ii*M0jj) in {q.fmt(6)}; raw {gm.fmt(5)}; mollification abs error <{float(deviation):.5f}')
    assert q.hi-q.lo<1000000*S  # sanity finite

# Exact support census: selected cell at (j-3/2)log3, central y in
# [-1/4-eps,1/4+eps], selected bump width ell=1/10000.
# The only valid signed prime powers are (j,n,sign)=(0,5,+),(1,2,+),
# (2,2,-),(3,5,-). For each n, all other centers are separated by >1/2.
# We certify the finite census with exact rational log enclosures.
def log_int(n):
    b=n.bit_length()-1;return b*LOG2+log_q(F(n,2**b))
ln3=log_int(3)
pps=[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]
shifts=[]
width=F(1,4)+eps+F(1,20000)
for j in range(4):
    c=ln3.scale(F(2*j-3,2))
    for n in pps:
        for sign in [-1,1]:
            z=c+log_int(n)*sign
            if z.lo>-width*S and z.hi<width*S:
                shifts.append((j,n,sign))
            else:
                assert z.lo>width*S or z.hi<-width*S,(j,n,sign,z.fmt())
assert shifts==[(0,5,1),(1,2,1),(2,2,-1),(3,5,-1)], shifts
assert (ln3/F(2)-Iv(F(1,4)+eps+F(1,20000))).lo>0
# w(t)<2 on the selected-central separation, certified using exp bounds.
# For t>=log3/2-1/4-eps-ell/2 > log(4/3)?? not quite;
# Instead use t>29/100 and w(t) decreasing, w(29/100)<2.
t0=ln3/F(2)-Iv(F(1,4)+eps+F(1,20000))
assert t0.lo>F(29,100)*S
# exp(-t/2)/(1-exp(-2t)) <2 at t=29/100 by elementary
# exp(-.145)<1/(1+.145) and exp(-.58)<1/(1+.58)
# crude inequalities insufficient; use alternating exp for t<=1.
eh=exp_neg_small(F(29,200));et=exp_neg_small(F(29,50))
ww=eh/(Iv(1)-et)
assert ww.hi<2*S,(ww.fmt(),t0.fmt())

# Conservative but genuine selected/central Weil cross estimate.
# For arbitrary y in each parity two-space, V proves ||y'||/||y||<38,46;
# Sobolev gives ||y||infty <= sqrt(2r)||y||, ||y||1<=sqrt(2a)||y||.
# Narrow selected-cell L1 <= sqrt(ell) * abs-coeff-sum /sqrt(Nsigma).
# Prime support census above implies exactly four selected-to-central
# translation overlaps: n=5 outside cells, n=2 inside cells.
# With s=4/sqrt3, N_±=(44∓16sqrt3)/3, coefficients 2s,2(s+2),
# and inner abs coeff |1∓s|, the exact interval-sup bounds below hold.
root3=sqrt_up(F(3));root3lo=sqrt_lo(F(3))
slo=F(4)/root3;shi=F(4)/root3lo
Npluslo=(F(44)-16*root3)/3
Nminuslo=(F(44)+16*root3lo)/3
k2=(log_int(2)/Iv(sqrt_lo(F(2)),sqrt_up(F(2)))).bounds()[1]
k5=(log_int(5)/Iv(sqrt_lo(F(5)),sqrt_up(F(5)))).bounds()[1]
assert k2<F(1,2) and k5<F(3,4)
for parity,r,desired in [(0,38,F(2,15)),(1,46,F(3,25))]:
    nlo=Npluslo if parity==0 else Nminuslo
    cinner=shi-1 if parity==0 else shi+1
    ctotal=2*shi if parity==0 else 2*(shi+2)
    # positive upper on normalized full source pairing:
    # prime max = sqrt(ell/n) (2 k5 + 2 cinner k2) sqrt(2r)
    # Gamma max = 2 sqrt(ell/n) ctotal sqrt(1/2+2 eps).
    ratio=F(1,100)/sqrt_lo(nlo)
    upper=ratio*((2*k5+2*cinner*k2)*sqrt_up(F(2*r))
                 +2*ctotal*sqrt_up(F(1,2)+2*eps))
    assert upper<desired,(parity,float(upper),desired)
    # mu_sigma>7, so compressed Schur correction <=cross²/7.
    schur=desired**2/7
    print('CERTIFIED',('even' if parity==0 else 'odd'),
          'ALL two-mode source cross |q(u,y)|/||y|| <',desired,
          'computed upper',str(float(upper)),
          '; compressed R operator <',schur)

# Rigorous Loewner two-by-two intervals on the genuinely smoothed matrix.
# A Hermitian 2x2 is positive if both diagonals and determinant are positive.
# Entry enclosures use independent worst-case intervals, so this is rigorous
# even though the intervals are correlated through a common mollifier.
def rigorous_2x2_spd(diag0,diag1,off):
    d0=diag0.bounds()[0];d1=diag1.bounds()[0]
    o=max(abs(off.bounds()[0]),abs(off.bounds()[1]))
    return d0>0 and d1>0 and d0*d1>o*o,(float(d0),float(d1),float(o),float(d0*d1-o*o))

for par,lo,hi,schur in [(0,F(17,20),F(31,20),F(4,1575)),
                         (1,F(23,20),F(9,5),F(9,4375))]:
    # Gamma - lo M is SPD: all Gamma Rayleigh values >lo.
    lower={ij:G_MOL[(par,*ij)]-lo*M_MOL[(par,*ij)]
           for ij in [(0,0),(0,1),(1,1)]}
    success,details=rigorous_2x2_spd(lower[(0,0)],lower[(1,1)],lower[(0,1)])
    assert success,(par,'lower gamma',details)
    # hi M - Gamma is SPD: all Gamma Rayleigh values <hi.
    upper={ij:hi*M_MOL[(par,*ij)]-G_MOL[(par,*ij)]
           for ij in [(0,0),(0,1),(1,1)]}
    success,details2=rigorous_2x2_spd(upper[(0,0)],upper[(1,1)],upper[(0,1)])
    assert success,(par,'upper gamma',details2)
    first=F(6)/(F(6)+hi)
    last=(F(6)+schur)/(F(6)+lo)
    assert first<F(1) and last<F(1)
    print('CERTIFIED',('even' if par==0 else 'odd'),
          'all two-dimensional normalized Gamma Rayleigh in (',lo,',',hi,')',
          '; Loewner determinant lower checks:',details,details2)
    print('CERTIFIED',('even' if par==0 else 'odd'),
          'B_F Ritz pencil has BOTH eigenvalues in (',first,',',last,')',
          'approx (',float(first),',',float(last),')')

print('PASS W: only original prime 2 and 5 enter selected-to-central couplings')
print('OPEN: interval upper spectral bound on the full complement, eta_N, RH')
