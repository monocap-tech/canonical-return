#!/usr/bin/env python3
"""CR-ALP-043: exact outward fixed-point band plus uniform tail budget.
Classical coefficient, digamma Euler-Maclaurin, and predecessor form theorems
are analytic premises, not Lean. No floating point in acceptance.
"""
from fractions import Fraction as F
from math import factorial,isqrt,comb
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
S=10**36; K=15000; CUTS=(172,242,300); a=F(11,10); sa=F(177,50)
def ceildiv(n,d):return -((-n)//d)
class I:
 __slots__=('l','h')
 def __init__(self,l,h=None,raw=False):
  if raw:self.l,self.h=l,h
  elif type(l) is int and h is None:self.l=self.h=l*S
  else:
   l=F(l);h=l if h is None else F(h);self.l=l.numerator*S//l.denominator;self.h=ceildiv(h.numerator*S,h.denominator)
  assert self.l<=self.h
 def __add__(x,y):
  y=iv(y);return I(x.l+y.l,x.h+y.h,True)
 __radd__=__add__
 def __neg__(x):return I(-x.h,-x.l,True)
 def __sub__(x,y):return x+-iv(y)
 def __rsub__(x,y):return iv(y)+-x
 def __mul__(x,y):
  if type(y) is int:return I(min(x.l*y,x.h*y),max(x.l*y,x.h*y),True)
  y=iv(y);z=(x.l*y.l,x.l*y.h,x.h*y.l,x.h*y.h);return I(min(z)//S,ceildiv(max(z),S),True)
 __rmul__=__mul__
 def __truediv__(x,y):
  if type(y) is int and y>0:return I(x.l//y,ceildiv(x.h,y),True)
  y=iv(y)
  if y.h<0:return (-x)/(-y)
  assert y.l>0,(y.l,y.h)
  z=[(u*S,v) for u in (x.l,x.h) for v in (y.l,y.h)]
  return I(min(n//d for n,d in z),max(ceildiv(n,d) for n,d in z),True)
 def __rtruediv__(x,y):return iv(y)/x
 def sq(x):
  lo=0 if x.l<=0<=x.h else min(x.l*x.l,x.h*x.h)
  return I(lo//S,ceildiv(max(x.l*x.l,x.h*x.h),S),True)
 def up(x):return F(max(abs(x.l),abs(x.h)),S)
 def sqrt(x):
  assert x.l>=0;return I(isqrt(x.l*S),isqrt(x.h*S)+1,True)
 def widen(x,n):return I(x.l-n,x.h+n,True)
 def mid(x):return F(x.l+x.h,2*S)
def iv(x):return x if isinstance(x,I) else I(x)
ZERO=I(0); ONE=I(1); QUARTER=I(F(1,4)); HALF=I(F(1,2))
# Truncation errors for reduced arguments. 48 terms each for log and atan.
assert 2*F(7,20)**97/(97*(1-F(7,20)**2))<F(1,S)
assert F(1,3)**97/97<F(1,S)
assert F(4)**64/factorial(64)<F(1,S)
def log_series(x):
 u=(x-1)/(x+1);assert u.up()<F(7,20);v=u.sq();t=u;z=ZERO
 for n in range(48):z+=t/(2*n+1);t=t*v
 return (2*z).widen(1)
LOG2=log_series(I(2))
def log(x):
 x=iv(x);assert x.l>0;k=0
 # range reduction by powers of 2, midpoint avoids boundary interval trouble
 while x.l+x.h>4*S:x=x/2;k+=1
 while x.l+x.h<2*S:x=x*2;k-=1
 return log_series(x)+k*LOG2
def atan(x):
 x=iv(x);assert x.up()<=F(1,3);v=-x.sq();t=x;z=ZERO
 for n in range(48):z+=t/(2*n+1);t=t*v
 return z.widen(1)
PI=16*atan(I(F(1,5)))-4*atan(I(F(1,239)))
assert PI.l>3*S and PI.h*7<22*S
TWO_PI=2*PI; NUF=PI/I(2*a)
def sincos(x):
 n=round(x.mid()/TWO_PI.mid());r=x-n*TWO_PI;assert r.up()<4
 q=-r.sq();s=r;c=ONE;ts=r;tc=ONE
 for m in range(1,32):
  ts=ts*q/((2*m)*(2*m+1));tc=tc*q/((2*m-1)*(2*m));s+=ts;c+=tc
 # cos through degree 62, sin through degree 63. Taylor error <= 4^64/64! for both.
 return s.widen(1),c.widen(1)
# Complex interval digamma by Euler-Maclaurin, order 12. Classical formula:
# psi(z)=log(z)-1/(2z)-sum B_(2m)/(2m z^(2m)) + remainder.
# sup |periodic B_24| <=4*24!/6^24; integral bound by real shift >=32,
# or by Im z >=120 and splitting t at Im z. Both <1e-31.
assert F(4*factorial(24),6**24*24*32**24)<F(1,10**31)
assert F(4*factorial(24),6**24*120**24)*(1+F(1,24))<F(1,10**31)
BER=[F(1)]
for n in range(1,25):BER.append(-sum(F(comb(n+1,k))*BER[k] for k in range(n))/(n+1))
assert BER[2]==F(1,6) and BER[4]==-F(1,30) and BER[24]==-F(236364091,2730)
def cmul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def inv(x,y):den=x.sq()+y.sq();return x/den,-y/den
def psi(x,y=ZERO):
 x,y=iv(x),iv(y);assert x.l>0 and y.l>=0
 n=0 if y.l>=120*S else 32;xx=x+n
 assert xx.l>=32*S or y.l>=120*S
 lr=log(xx.sq()+y.sq())/2
 li=atan(y/xx) if n else PI/2-atan(xx/y)
 z=inv(xx,y);u=cmul(z,z);p=u
 re,im=lr-z[0]/2,li-z[1]/2
 for m in range(1,13):
  b=I(BER[2*m]/(2*m));re-=b*p[0];im-=b*p[1];p=cmul(p,u)
 for m in range(n):
  h=inv(x+m,y);re-=h[0];im-=h[1]
 return re.widen(10**5),im.widen(10**5) # 1e-31
# Pinned source intervals and saved rational basis. Recheck with predecessor separately.
HF=ROOT/'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv';IF=ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt'
assert hashlib.sha256(HF.read_bytes()).hexdigest()=='7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c'
assert hashlib.sha256(IF.read_bytes()).hexdigest()=='6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9'
P={}
for ln in IF.read_text().splitlines():
 k,l,h=ln.split();P[k]=I(F(int(l),10**25),F(int(h),10**25))
prior=json.loads((ROOT/'CR-ALP-041_energy_conditioned_result.json').read_text())
N0=I(a).sqrt(); q=P['EINV'];assert q.up()<F(2,5)
# 32 exponential terms; beta>=1/2 and tail m>=32. All the omitted
# exponential correction sums below have term factor <=10; hence radius <1e-36.
assert 10*F(2,5)**129/(1-F(2,5)**4)<F(1,S)
BETA=[I(F(4*m+1,2)) for m in range(32)]; B2=[b.sq() for b in BETA]; EX=[q]
q4=q.sq().sq()
for m in range(1,32):EX.append(EX[-1]*q4)
psi0=psi(QUARTER)[0]; S1=1-LOG2; S2=2*psi(I(F(5,4)))[0]-psi(ONE)[0]-psi(I(F(3,2)))[0]
NS=(2,3,4,5,7,8,9); primes=(2,3,2,5,7,2,3)
TS=[log(I(n)) for n in NS]; ROOTS=[I(n).sqrt() for n in NS]; W=[log(I(p))/r for p,r in zip(primes,ROOTS)]
assert sum(w.up() for w in W)<sa
NUS={i:i*NUF for i in range(1,15)}
# Exact small-index sine-centering constants: R_s(v)=sum(1-sigma exp(-2a beta))/(beta²+v²).
def Rbase(v):return psi(QUARTER,v/2)[1]/(2*v)
def Rweighted(v,sgn):return Rbase(v)-sgn*sum(e/(b+v.sq()) for e,b in zip(EX,B2)).widen(1)
AI={}
for i,u in NUS.items():
 theta=sum(w*sincos(u*t)[0] for w,t in zip(W,TS));AI[i]=u*Rweighted(u,(-1)**i)+theta
# Tail envelope for the *specific* pole-centered coefficients, not just projected norms.
def remainder_bound(ids,co,y,k,N):
 # CR-ALP-038 formulas apply also to an ordinary sine component: m=0 is
 # removed by the chosen pole subtraction (beta_0²-1/4=0).
 D=sum(co[i]*NUS[i] for i in ids).up()
 C={s:sum(co[i].up()*NUS[i].up()**s for i in ids) for s in (0,2,3)}
 A=(2*sa/a+F(201,200)*F(22,7)/(2*a))*(2*a/3)**2*D
 B=(sa*C[0]/(2*a)+8*sa*C[2]/a+F(201,200)*C[3]/a)*(2*a/3)**3
 ri=A/isqrt(3*K**3)+B/isqrt(5*K**5)
 B0,B1=(F(6,5),F(3,10)) if k==1 else (F(3,5),F(3,5))
 # CR-ALP-039 centered log envelope; monotone log²(j/N)/j² for j>K>eN.
 assert K>3*N
 l=log(I(F(K,N))).up(); l2=I((l*l+2*l+2)/K).sqrt().up()
 C2=F(121,225)*(F(12001,4000)*B0+F(11,7)*B1);C3=8*a**3*B0/945
 rg=F(7,10)*B0*l2+C2/isqrt(3*K**3)+C3/isqrt(5*K**5)
 # Exact b-translation centered formula gives this residual bound directly.
 rp=2*B0*sa*(2*a/3)/isqrt(K)+2*B1*sa*(2*a/3)**2/isqrt(3*K**3)
 return ri+(y*P[f'gamma_{k}']).up()*(rg+rp)

def run():
 results=[]
 for bl in prior['blocks']:
  ids=bl['indices'];k=ids[-1];sgn=(-1)**k;s=-sgn
  b0=(P['EHALF']+P['EHI'])/2 if k==1 else -(P['EHALF']-P['EHI'])/2
  b1=-(P['EHALF']-P['EHI'])/4 if k==1 else (P['EHALF']+P['EHI'])/4
  ep=(b0+2*b1)/2;em=(b0-2*b1)/2
  E0=ep*(2*a)+em*(1-P['EINV'].sq())-2*b0
  # E_m = [b0/4+beta b1]/[beta(beta²-.25)] - s exp*(beta b0-b1)/(beta²-.25), m>=1.
  EC=[s*e*(b*b0-b1)/(bb-QUARTER) for b,bb,e in zip(BETA[1:],B2[1:],EX[1:])]
  Einf=E0+b1*S1+b0/4*S2-sum(EC).widen(1)
  t0=(log(NUF*CUTS[0]/2)-psi0)/2-Einf/b0
  shift={N:log(I(F(N,CUTS[0])))/2 for N in CUTS};shift[CUTS[0]]=ZERO
  coords=[];d=[];ys=[]
  for v,bu in zip(bl['basis_rows'],bl['budgets']):
   x=[I(F(z)) for z in v];y=x[-1];co={i:z for i,z in zip(ids,x)};co[k]-=sum(co[i]*P[f'lam_{i}'] for i in ids[:-1]);coords.append(co);ys.append(y);d.append(I(F(bu['d'])))
  # Precompute transformed sine constants directly; preserve signed cancellation.
  trace={N:I(0) for N in CUTS}; max_width=0; nrows=0
  trigstep=[sincos(2*NUF*t) for t in TS]
  trig=None
  for j in range(CUTS[0]+1+(k==2),K+1,2):
   v=j*NUF;v2=v.sq();den=v2+QUARTER
   psir,psii=psi(QUARTER,v/2);rv=psii/(2*v);Lam=psir-psi0
   epart=sum(e/(bb+v2) for e,bb in zip(EX,B2)).widen(1)
   Rj=rv-sgn*epart
   if nrows%16==0:trig=[sincos(v*t) for t in TS]
   else:trig=[(u*c+v0*s,v0*c-u*s) for (u,v0),(s,c) in zip(trig,trigstep)]
   ss=sum(w*t[0] for w,t in zip(W,trig));cc=sum(w*t[1] for w,t in zip(W,trig))
   # Actual source-signed sine coefficient after subtracting t_i c_j.
   CS={i:2/I(a)*(v*AI[i]/den-(v*AI[i]-NUS[i]*(v*Rj+ss))/(v2-NUS[i].sq())) for i in ids}
   cj=2*b0*v/(N0*den)
   # Actual b Gamma coefficient with the convergent E_m series summed by partial fractions.
   R1=rv-1/den;S3=Lam/(2*v2)-2/den
   Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v2) for e,bb in zip(EC,B2[1:])).widen(1)
   Cb=Lam*(cj-b0/(N0*v))-2*v/N0*Es-t0*cj-2/(N0*den)*(v*b0*cc-b1*ss)
   for co,y,di in zip(coords,ys,d):
    yg=y*P[f'gamma_{k}'];base=sum(co[i]*CS[i] for i in ids)-yg*Cb
    for N in CUTS:
     if j>N:
      z=base+yg*shift[N]*cj;trace[N]+=z.sq()/di
      max_width=max(max_width,z.h-z.l)
   nrows+=1
   if nrows%1000==0:print('parity',k,'processed',nrows,flush=True)
  for N in CUTS:
   tr=F(trace[N].h,S);tail=sum(remainder_bound(ids,co,y,k,N)**2/F(bu['d']) for co,y,bu in zip(coords,ys,bl['budgets']))
   # Round outward to millionths for compact replayable exact acceptance.
   def rup(x):return F(ceildiv(x.numerator*10**6,x.denominator),10**6)
   br=rup(tr);ta=rup(tail);rho=br+ta;mu={172:F(1,4),242:F(3,5),300:F(3,4)}[N]
   print('parity',k,'cutoff',N,'band upper',float(br),'remainder upper',float(ta),'rho upper',float(rho),'mu',mu,flush=True)
   results.append(dict(parity=k,N=N,band_end=K,band_rows=nrows,band_upper=str(br),tail_upper=str(ta),rho_upper=str(rho),mu=str(mu),pass_positive=rho<mu,max_coefficient_interval_width_units=str(max_width)))
 # Explicit acceptance only if all parity comparisons pass for these cutoffs.
 for N,target in [(242,F(2,5)),(300,F(2,5))]:
  assert all(F(r['rho_upper'])<target for r in results if r['N']==N),(N,results)
  assert target<({242:F(3,5),300:F(3,4)}[N])
 out=dict(stage='CR-ALP-043',status='Exact outward-integer coefficient band and analytic infinite remainder PASS at 242,300, conditional on continuum premises',results=results,not_claimed='E172 mixed, E242 or E300 full finite section, full collar, arbitrary support, Lean, RH')
 (ROOT/'CR-ALP-043_exact_band_result.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS exact band+remainder inequalities; no full collar or RH claim',flush=True)
if __name__=='__main__':run()
