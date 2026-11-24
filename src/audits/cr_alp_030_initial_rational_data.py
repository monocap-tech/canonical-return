#!/usr/bin/env python3
"""CR-ALP-030: exact rational INPUT GENERATOR for original signed eight-mode certificate.

This extends the CR-ALP-029 six interior w3..w8 functions by the true
physical pole-null boundary anchors u1,u2. Retains the EXACT original
seven negative paired prime powers, digamma gamma, and physical metric.

Mathematical inputs (not Lean-proved here): source EF, digamma translation-
difference PSD expansion, exact sine and exponential Laplace identities,
Euler-Maclaurin Euler constant bound, and exact continuous Fourier form.
This program produces a rigorously outward-contained 10^-25 input file for
an independently compiled checked C++ interval summation and final LDL.
It audits physical mass and complete signed prime entries. It does NOT
certify all 15k modes, arbitrary collar, RH, or a Lean theorem.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import importlib.util,sys,json

HERE=Path(__file__).resolve().parent
PREV=HERE/'cr_alp_028_prime_interval_and_gram_audit.py'
if not PREV.exists(): PREV=HERE.parent/'cr_alp_028'/'cr_alp_028_prime_interval_and_gram_audit.py'
spec=importlib.util.spec_from_file_location('original_prime',PREV)
prime=importlib.util.module_from_spec(spec);spec.loader.exec_module(prime)
I=prime.I;A0=prime.A;PI=prime.pi

def log_interval_pos(z, terms=65):
 z=I(z);assert 1<=z.lo<=z.hi<=2
 u=(z-I(1))/(z+I(1)); lo=I(0);powu=u
 for j in range(terms):lo=lo+2*powu/I(2*j+1);powu=powu*u*u
 tail=2*powu/(I(2*terms+1)*(I(1)-u*u))
 return I(lo.lo,(lo+tail).hi)

def exp_rational_pos(z,terms=75):
 z=F(z);assert F(0)<=z<2
 s=sum((z**j/F(factorial(j)) for j in range(terms+1)),F(0))
 tail=(z**(terms+1)/factorial(terms+1))/(1-z/F(terms+2))
 return I(s,s+tail)

n=256
Hn=sum((F(1,k) for k in range(1,n+1)),F(0))
lnn=prime.log_int(n)
D4=F(1,2*n)-F(1,12*n**2)+F(1,120*n**4)
D6=D4-F(1,252*n**6)
EULER=I((I(Hn)-lnn-I(D4)).lo,(I(Hn)-lnn-I(D6)).hi)
ln_pi=prime.log_int(2)+log_interval_pos(PI/I(2))
g0=-EULER-PI/I(2)-I(3)*prime.log_int(2)-ln_pi
assert g0.lo>F(-6) and g0.hi<F(-5)

SCALE=10**42
class FI:
 __slots__=('lo','hi')
 @staticmethod
 def make(lo,hi):
  v=object.__new__(FI);v.lo,v.hi=lo,hi;return v
 def __init__(self,x=0):
  if isinstance(x,FI): self.lo,self.hi=x.lo,x.hi;return
  if isinstance(x,I):
   self.lo=x.lo.numerator*SCALE//x.lo.denominator
   self.hi=-((-x.hi.numerator*SCALE)//x.hi.denominator)
   return
  x=F(x);self.lo=x.numerator*SCALE//x.denominator
  self.hi=-((-x.numerator*SCALE)//x.denominator)
 def __add__(self,b):
  b=FI(b);return FI.make(self.lo+b.lo,self.hi+b.hi)
 __radd__=__add__
 def __neg__(self):return FI.make(-self.hi,-self.lo)
 def __sub__(self,b):return self+-FI(b)
 def __rsub__(self,b):return FI(b)+-self
 def __mul__(self,b):
  b=FI(b); vals=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
  return FI.make(min(vals)//SCALE,-((-max(vals))//SCALE))
 __rmul__=__mul__
 def __truediv__(self,b):
  b=FI(b);assert b.lo>0,'denominator nonpositive'
  a=[n*SCALE//d for n in (self.lo,self.hi) for d in (b.lo,b.hi)]
  c=[-((-n*SCALE)//d) for n in (self.lo,self.hi) for d in (b.lo,b.hi)]
  return FI.make(min(a),max(c))
 def __repr__(self):return f'{self.lo}/{SCALE} .. {self.hi}/{SCALE}'

A=FI(A0);# Explicit rigorous sqrt(11/10) interval by integer isqrt
from math import isqrt
rt_scale=10**45
sq=isqrt(11*rt_scale**2//10) # floor sqrt(11/10)*rt_scale
sqrtA=FI(I(F(sq,rt_scale),F(sq+1,rt_scale)))

EHALF=FI(exp_rational_pos(F(11,20))); E=EHALF*EHALF
EINV=FI(1)/E; EHI=FI(1)/EHALF
q=EINV  # e^-a, so exp(-2*a*beta)=q^(4k+1)
assert q.lo>3*SCALE//10 and q.hi<2*SCALE//5
SinhA=(E-EINV)/FI(2)
CoshHalf=(EHALF+EHI)/FI(2)
SinhHalf=(EHALF-EHI)/FI(2)

NU={j:FI(prime.nu[j]) for j in range(1,9)}
lam={j:FI(prime.lam(j)) for j in range(3,9)}
B={1:A+SinhA,2:SinhA-A}
c={1:FI(2)*NU[1]*CoshHalf/(sqrtA*(NU[1]*NU[1]+FI(F(1,4)))),
   2:-FI(2)*NU[2]*SinhHalf/(sqrtA*(NU[2]*NU[2]+FI(F(1,4))))}
gamma={k:c[k]/B[k] for k in (1,2)}
eta={k:c[k]*c[k]/B[k] for k in (1,2)}
assert all(B[k].lo>0 for k in (1,2))
assert all((FI(1)-eta[k]).lo>0 for k in (1,2))

# Source-authentic 7 weights and t=log(n), already outward certified.
weights={n:FI(prime.weights[n]) for n in prime.POWERS}
times={n:FI(prime.t[n]) for n in prime.POWERS}
sqrtn={n:FI(prime.sqrt_int(n)) for n in prime.POWERS}

# Construct physical basis by parity: three inner w_j plus u_k.
BLOCKS=[((3,5,7),1,'even physical'),((4,6,8),2,'odd physical')]
M={};H={}
for js,k,_ in BLOCKS:
 for i in js:
  for j in js:M[(i,j)]=FI(int(i==j))+lam[i]*lam[j]
  M[(i,k)]=M[(k,i)]=-lam[i]
 M[(k,k)]=FI(1)-eta[k]

# Odd physical k=2 uses sinh basis; even physical k=1 cosh.
def exp_pt(p,n):return sqrtn[n] if p==F(1,2) else FI(1)/sqrtn[n]
def exp_pa(p):return EHALF if p==F(1,2) else EHI

def Uexp_sine(j,p,n):
 """U(e_j, exp(p*x)) for t=log(n), EXACT integral of translations."""
 v=NU[j];tr=prime.trig(prime.nu[j]*prime.t[n],'sin')
 cr=prime.trig(prime.nu[j]*prime.t[n],'cos')
 ep=exp_pa(p); em=FI(1)/ep; tau=exp_pt(p,n)
 sgn=(-1)**j
 term1=em*(v*FI(cr)-FI(p)*FI(tr)) - FI(sgn)*v*ep/tau + FI(sgn)*ep*(-FI(p)*FI(tr)-v*FI(cr)) +v*tau*em
 return term1/(sqrtA*(v*v+FI(F(1,4))))

def Uexp_exp(p,r,n):
 """U(exp(px),exp(rx)); translate with t=log n, fixed support."""
 t=times[n];taup=exp_pt(p,n);taur=exp_pt(r,n)
 if p+r==0:return (FI(2)*A-t)*(FI(1)/taup+FI(1)/taur)
 # exact interval for exponent (p+r)=+/-1, exp[(p+r)a], exp[-(p+r)a]
 Ep=E if p+r==1 else EINV
 Em=FI(1)/Ep
 # exp[(p+r)*t] = n^(p+r), n or 1/n
 shift=FI(n) if p+r==1 else FI(F(1,n))
 overlap=((Ep-Em*shift)/FI(p+r)) if p+r>0 else ((Em*shift-Ep)/FI(-(p+r)))
 return (FI(1)/taup+FI(1)/taur)*overlap

def U_sine_b(j,k,n):
 plus=Uexp_sine(j,F(1,2),n);minus=Uexp_sine(j,F(-1,2),n)
 return (plus+(minus if k==1 else -minus))/FI(2)

def U_b_b(k,n):
 pp=Uexp_exp(F(1,2),F(1,2),n)
 mm=Uexp_exp(F(-1,2),F(-1,2),n)
 pm=Uexp_exp(F(1,2),F(-1,2),n)
 return (pp+mm+(FI(2)*pm if k==1 else -FI(2)*pm))/FI(4)

prime_Psine={(i,j):FI(prime.P(i,j)) for i in range(1,9) for j in range(1,9) if i%2==j%2}

def P_sine_b(j,k):
 return -sum((weights[n]*U_sine_b(j,k,n) for n in prime.POWERS),FI(0))
def P_bb(k):
 return -sum((weights[n]*U_b_b(k,n) for n in prime.POWERS),FI(0))

for js,k,label in BLOCKS:
 pb={j:P_sine_b(j,k) for j in (k,*js)}
 pbb=P_bb(k)
 for i in js:
  for j in js:
   # source previous six-mode audited prime intervals
   H[(i,j)]=FI(g0)*M[(i,j)]+FI(prime.Pw(i,j))
  H[(i,k)]=H[(k,i)]=FI(g0)*M[(i,k)] + prime_Psine[(i,k)]-lam[i]*prime_Psine[(k,k)] -gamma[k]*(pb[i]-lam[i]*pb[k])
 H[(k,k)]=FI(g0)*M[(k,k)] +prime_Psine[(k,k)]-FI(2)*gamma[k]*pb[k]+gamma[k]*gamma[k]*pbb

# Gamma beta matrix on sines and b_k. J=2M/beta-K_beta (K_beta has e^-beta|x-y| kernel).
def Iexp(z):
 if z==0:return FI(2)*A
 assert z in (F(1),F(-1))
 return FI(2)*SinhA # I(exp(+-x)) = 2sinh(a)

def Iprime(z):
 if z==0:return FI(0)
 assert z in (F(1),F(-1))
 # derivative I'(z)=2a cosh(az)/z -2sinh(az)/z²
 cosh=(E+EINV)/FI(2)
 v=FI(2)*A*cosh-FI(2)*SinhA
 return v if z==1 else -v

def Lbeta(p,r,beta):
 """Double half-plane Laplace ∫_{x>y}e^{px+ry-beta(x-y)}dy dx."""
 z=r+beta
 if z==0:return Iprime(p-beta)+A*Iexp(p-beta)
 # At beta=1/2, z=r+beta ∈{0,1}; when z=1 need exp(-a)
 if z==1:ep=EINV
 return (Iexp(p+r)-ep*Iexp(p-beta))/FI(z)

def Kbb_beta(k,beta,qexp):
 ps=(F(1,2),F(-1,2))
 if beta==F(1,2):
  val={(p,r):Lbeta(p,r,beta)+Lbeta(r,p,beta) for p in ps for r in ps}
 else:
  val={}
  # coefficient qexp=e^-2a beta; for 1-exp[-2a(beta-r)] = 1-qexp*e^(2ar), e^(2ar)=E or 1/E
  val={}
  for p in ps:
   for r in ps:
    ep=E if p+r==1 else EINV if p+r==-1 else FI(1)
    em=FI(1)/ep
    er=E if r==F(1,2) else EINV
    term1=FI(2)*FI(beta)*Iexp(p+r)/FI(beta*beta-p*p)
    term2=em*(FI(1)-qexp*er)/FI((beta+p)*(beta-r))
    # term3 exp(-2a(beta+r)) = qexp / exp(2ar)
    term3=ep*(FI(1)-qexp/er)/FI((beta-p)*(beta+r))
    val[p,r]=term1-term2-term3
 if k==1:
  return (val[F(1,2),F(1,2)]+val[F(-1,2),F(-1,2)]+val[F(1,2),F(-1,2)]+val[F(-1,2),F(1,2)])/FI(4)
 return (val[F(1,2),F(1,2)]+val[F(-1,2),F(-1,2)]-val[F(1,2),F(-1,2)]-val[F(-1,2),F(1,2)])/FI(4)


def K_sine_b(j,k,beta,qexp):
 v=NU[j];db=FI(beta*beta)+v*v
 m={}
 for p in (F(1,2),F(-1,2)):
  ep=exp_pa(p);em=FI(1)/ep
  mom=v*(em-FI((-1)**j)*ep)/(sqrtA*(v*v+FI(F(1,4))))
  # A_left=e^-pa*(1-exp[-2a(beta-p)])/(beta-p)
  # A_right=e^pa*(1-exp[-2a(beta+p)])/(beta+p)
  fac=E if p==F(1,2) else EINV # exp(2ap)
  ld=FI(2)*A*em if beta==p else em*(FI(1)-qexp*fac)/FI(beta-p)
  rd=FI(2)*A*ep if beta==-p else ep*(FI(1)-qexp/fac)/FI(beta+p)
  m[p]=(FI(2)*FI(beta)*mom+v/sqrtA*(ld-FI((-1)**j)*rd))/db
 return (m[F(1,2)]+(m[F(-1,2)] if k==1 else -m[F(-1,2)]))/FI(2)

def K_sine(i,j,beta,qexp):
 ni=NU[i];nj=NU[j]
 if i%2!=j%2:return FI(0)
 di=FI(beta*beta)+ni*ni; dj=FI(beta*beta)+nj*nj
 ans=FI(2)*FI(beta)/di if i==j else FI(0)
 ans+=FI(2)*ni*nj*(FI(1)-FI((-1)**i)*qexp)/(A*di*dj)
 return ans


# Outward convert 10^-42 integer intervals to checked 10^-25 integer intervals
# WITHOUT losing containment; physical source constants remain enclosed.
INITIAL_SCALE = 10**25
OUTWARD_FACTOR = SCALE // INITIAL_SCALE
assert SCALE == 10**42 and OUTWARD_FACTOR == 10**17
values={k:v for k,v in [('A',A),('SQRT_A',sqrtA),('E',E),
  ('EINV',EINV),('EHALF',EHALF),('EHI',EHI),('g0',FI(g0)),
  ('q',q),('SinhA',SinhA)]}
for j in range(1,9):values[f'nu_{j}']=NU[j]
for j in range(3,9):values[f'lam_{j}']=lam[j]
for k in (1,2):values[f'gamma_{k}']=gamma[k]
for (i,j),v in M.items():values[f'M_{i}_{j}']=v
for (i,j),v in H.items():values[f'H_{i}_{j}']=v
out=HERE/'CR-ALP-030_initial_rational_intervals_1e25.txt'
with out.open('w') as f:
 for key,v in sorted(values.items()):
  lo=v.lo//OUTWARD_FACTOR
  hi=-((-v.hi)//OUTWARD_FACTOR)
  assert lo<=hi
  f.write(f'{key} {lo} {hi}\n')
print('CR-ALP-030 exact source-initial intervals: PASS',len(values),'constants')
print('Saved:',out)
