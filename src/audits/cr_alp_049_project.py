#!/usr/bin/env python3
"""CR-ALP-049: rational trial inverse and original signed residual intervals.
Classical coefficient and Euler-Maclaurin identities are analytical inputs.
Exact frozen witness replay without midpoint optimization. No Lean proof.
"""
from fractions import Fraction as F
from math import factorial,isqrt,comb
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"cr_alp_045"))
from cr_alp_045_layered_gamma_audit import lower
mu,_=lower(85)
gamma=F(72,25)
assert mu+F(21,8)>gamma
print("gamma floor",gamma,"scalar remainder",gamma-F(21,8),flush=True)
S=10**100; a=F(11,10)
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
# Truncation errors for reduced arguments. 400 terms each for log and atan.
assert 2*F(7,20)**801/(801*(1-F(7,20)**2))<F(1,S)
assert F(1,3)**801/801<F(1,S)
assert F(4)**300/factorial(300)<F(1,S)
def log_series(x):
 u=(x-1)/(x+1);assert u.up()<F(7,20);v=u.sq();t=u;z=ZERO
 for n in range(400):z+=t/(2*n+1);t=t*v
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
 for n in range(400):z+=t/(2*n+1);t=t*v
 return z.widen(1)
PI=16*atan(I(F(1,5)))-4*atan(I(F(1,239)))
assert PI.l>3*S and PI.h*7<22*S
TWO_PI=2*PI; NUF=PI/I(2*a)
def sincos(x):
 n=round(x.mid()/TWO_PI.mid());r=x-n*TWO_PI;assert r.up()<4
 q=-r.sq();s=r;c=ONE;ts=r;tc=ONE
 for m in range(1,150):
  ts=ts*q/((2*m)*(2*m+1));tc=tc*q/((2*m-1)*(2*m));s+=ts;c+=tc
 # cos through degree 798, sin through degree 799. Error <= 4^800/800!.
 return s.widen(1),c.widen(1)
# Euler-Maclaurin with derivative, shift 512 and 90 Bernoulli terms.
M=24;SHIFT=512;N=410;KTRIAL=205
assert F(4*factorial(2*M),6**(2*M)*2*M*SHIFT**(2*M))<F(1,10**92)
assert F(4*factorial(2*M),6**(2*M)*SHIFT**(2*M+1))<F(1,10**92)
BER=[F(1)]
for n in range(1,2*M+1):BER.append(-sum(F(comb(n+1,k))*BER[k] for k in range(n))/(n+1))
def cmul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def inv(x,y):den=x.sq()+y.sq();return x/den,-y/den
def psi_pair(x,y=ZERO):
 x,y=iv(x),iv(y);xx=x+SHIFT;assert y.up()/F(xx.l,S)<F(3,5)
 lr=log(xx.sq()+y.sq())/2;li=(atan(y/xx) if (y/xx).up()<=F(1,3) else atan(I(F(1,4)))+atan((y/xx-F(1,4))/(1+y/xx/4)));z=inv(xx,y);u=cmul(z,z);p=u
 re,im=lr-z[0]/2,li-z[1]/2;tr=z[0]+u[0]/2
 for m in range(1,M+1):
  b=I(BER[2*m]/(2*m));re-=b*p[0];im-=b*p[1];tr+=I(BER[2*m])*cmul(p,z)[0];p=cmul(p,u)
 for m in range(SHIFT):
  h=inv(x+m,y);re-=h[0];im-=h[1];tr+=cmul(h,h)[0]
 return re.widen(10**8),im.widen(10**8),tr.widen(10**8)
def exp_pos(x):
 assert 0<=x.l and x.up()<3
 t=ONE;z=ONE
 for n in range(1,400):t=t*x/n;z+=t
 assert F(3)**400/factorial(400)/(1-F(3,401))<F(1,S)
 return z.widen(1)
A=I(a);L=2*A;rt=A.sqrt();E=exp_pos(A);eh=exp_pos(A/2);em=1/eh;q=1/E;qt=QUARTER
assert q.up()<F(2,5)
psi0=psi_pair(qt)[0];g0=psi0-log(PI)
S1=1-LOG2;S2=2*psi_pair(I(F(5,4)))[0]-psi_pair(ONE)[0]-psi_pair(I(F(3,2)))[0];S4=PI.sq()/6-1-2*S1
K=256;assert 100*F(2,5)**(4*K+1)/(1-F(2,5)**4)<F(1,S)
B=[I(F(4*m+1,2)) for m in range(K)];B2=[b.sq() for b in B];EX=[q]
for m in range(1,K):EX.append(EX[-1]*q.sq().sq())
NS=(2,3,4,5,7,8,9);PS=(2,3,2,5,7,2,3);TS=[log(I(n)) for n in NS];W=[log(I(p))/I(n).sqrt() for p,n in zip(PS,NS)]
bankpath=ROOT/'CR-ALP-046_full_E85_intervals.json'
if not bankpath.exists():bankpath=ROOT.parent/'cr_alp_046/CR-ALP-046_full_E85_intervals.json'
raw=bankpath.read_bytes(); assert hashlib.sha256(raw).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
bank=json.loads(raw); oldS=bank['scale']; div=oldS//S; assert oldS==div*S
nu={j:j*NUF for j in range(1,N+1)}; trig={}; R={}; G={}; AI={}; R2={}
for j,v in nu.items():
 re,im,tre=psi_pair(qt,v/2);R[j]=im/(2*v);G[j]=re-log(PI);trig[j]=[sincos(v*t) for t in TS]; AI[j]=v*(R[j]-(-1)**j*sum(e/(b+v.sq()) for e,b in zip(EX,B2)).widen(1))+sum(w*t[0] for w,t in zip(W,trig[j]))
 R2[j]=im/(4*v*v*v)-tre/(8*v.sq())-(-1)**j*sum(e/(b+v.sq()).sq() for e,b in zip(EX,B2)).widen(1)
 if j%20==0:print('interval coefficients',j,flush=True)
def hss(i,j):
 u,v=nu[i],nu[j]
 if i==j:return G[i]-2*u.sq()/A*R2[i]-sum(w*((L-t)/A*c+s/(A*u)) for w,t,(s,c) in zip(W,TS,trig[i]))
 return -2/A*(v*AI[i]-u*AI[j])/(v.sq()-u.sq())
def certify(X):
 # integer midpoint minus maximum interval row radius; outward compression to 65 digits
 d=S//10**65;n=len(X);C=[[0]*n for _ in range(n)];r=[0]*n
 for i in range(n):
  for j in range(n):
   l=X[i][j].l//d;h=ceildiv(X[i][j].h,d);c=(l+h)//2;C[i][j]=c;r[i]+=max(c-l,h-c)
 for i in range(n):
  for j in range(i):assert C[i][j]==C[j][i]
 for i in range(n):C[i][i]-=max(r)
 prev=1;piv=[]
 for k in range(n):
  pivot=C[k][k];piv.append(pivot>0)
  if pivot<=0:break
  for i in range(k+1,n):
   for j in range(i,n):
    v,rem=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert rem==0;C[i][j]=C[j][i]=v
  prev=pivot
 return dict(dimension=n,positive_pivots=sum(piv),all_positive=len(piv)==n and all(piv),row_radius_units=str(max(r)),scale='1e65')
from functools import lru_cache
primeAI={j:sum(w*t[0] for w,t in zip(W,trig[j])) for j in nu}
@lru_cache(None)
def pss(i,j):
 u,v=nu[i],nu[j]
 if i==j:return -sum(w*((L-t)/A*ct+st/(A*u)) for w,t,(st,ct) in zip(W,TS,trig[i]))
 return -2/A*(v*primeAI[i]-u*primeAI[j])/(v.sq()-u.sq())
import mpmath as mp
import numpy as np
mp.mp.dps=85

def mpiv(x):return mp.mpf(x.l+x.h)/(2*S)
def matr(A):return mp.matrix([[mpiv(x) for x in row] for row in A])
def rat(A,den=10**55):return np.array([[I(F(int(mp.nint(A[i,j]*den)),den)) for j in range(A.cols)] for i in range(A.rows)],dtype=object)
def mat(A):return np.array(A,dtype=object)
def upper(A):return max(sum(x.up() for x in row) for row in A)
def compact(x):return [str(x.l),str(x.h)]
for bl in bank['blocks']:
 k=bl['parity'];ids=bl['ids'];n=len(ids);s=-(-1)**k
 H=mat([[I(l//div,ceildiv(h,div),True) for l,h in row] for row in bl['H']])
 b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4
 BM=(E-q)/2+(A if k==1 else -A);ep=(b0+2*b1)/2;en=(b0-2*b1)/2
 E0=ep*L+en*(1-q.sq())-2*b0;EC=[s*e*(b*b0-b1)/(bb-qt) for b,bb,e in zip(B[1:],B2[1:],EX[1:])]
 def c(j):v=nu[j];return 2*b0*v/(rt*(v.sq()+qt))
 def hb(j):
  v=nu[j];den=v.sq()+qt;Lam=G[j]-g0;R1=R[j]-1/den;S3=Lam/(2*v.sq())-2/den
  Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v.sq()) for e,bb in zip(EC,B2[1:])).widen(1)
  val=g0*c(j)+Lam*(c(j)-b0/(rt*v))-2*v/rt*Es
  for nn,w,(st,ct) in zip(NS,W,trig[j]):
   sn=I(nn).sqrt();bt=ep*sn+en/sn;val-=2*w/(rt*den)*(v*(bt+b0*ct)-b1*st)
  return val
 band=list(range(87 if k==1 else 86,KTRIAL+1,2));last=band[-1];js0=band[:-1];nb=len(js0)
 rows=mat([[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-c(k)/BM*hb(j) for i in ids] for j in band])
 V=rat(mp.cholesky(matr(H)).T**-1)
 EV=V.T@H@V
 err=upper(EV-mat([[I(int(i==j)) for j in range(n)] for i in range(n)]));assert err<F(1,10**12)
 print('parity',k,'energy congruence error',float(err),flush=True)
 U=rows@V
 Fmat=mat([[I(int(i==j)) if i!=nb else -c(js0[j])/c(last) for j in range(nb)] for i in range(nb+1)])
 Lraw=mat([[(I(gamma) if i==j else ZERO)+pss(i,j) for j in band] for i in band])
 LF=Fmat.T@Lraw@Fmat; CF=Fmat.T@U
 Y=rat(mp.inverse(matr(LF))@matr(CF))
 X=Fmat@Y
 Kmat=EV-(CF.T@Y+Y.T@CF-Y.T@LF@Y)
 Rraw=U-Lraw@X
 poles=mat([[c(j)] for j in band]);fit=rat(matr(poles.T@Rraw)/mpiv((poles.T@poles)[0,0]));res=Rraw-poles@fit
 for i in range(n):
  for j in range(i):Kmat[i,j]=Kmat[j,i]=I(max(Kmat[i,j].l,Kmat[j,i].l),min(Kmat[i,j].h,Kmat[j,i].h),True)
 # Certify a quantitative energy margin 1/4, 1/8, etc. by exact interval determinant test.
 cert=None
 for m in (F(1,2),F(1,4),F(1,8),F(1,16),F(1,32),F(0)):
  test=certify([[Kmat[i,j]-(I(m) if i==j else ZERO) for j in range(n)] for i in range(n)])
  if test['all_positive']:cert=test;break
 assert cert is not None
 rho=upper(res.T@res);assert rho<F(1,10**12)
 print('parity',k,'finite corrected energy margin',m,'pivots',cert['positive_pivots'],'in-band fitted residual Gram upper',float(rho),flush=True)
 # Explicit omitted band, same fitted pole subtraction globally.
 outer=list(range(207 if k==1 else 206,N+1,2))
 rows2=mat([[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-c(k)/BM*hb(j) for i in ids] for j in outer])
 U2=rows2@V
 R2band=U2-mat([[pss(i,j) for j in band] for i in outer])@X-mat([[c(j)] for j in outer])@fit
 # Exact orthogonal projection onto F_{85,410}: a lower bound on the full
 # residual Gram, independent of the arbitrary pole fit.
 rr=np.concatenate([res,R2band],axis=0);pp=mat([[c(j)] for j in band+outer]);pc=pp.T@rr;psq=(pp.T@pp)[0,0]
 GB=rr.T@rr-np.array([[pc[0,i]*pc[0,j]/psq for j in range(n)] for i in range(n)],dtype=object)
 KK=Kmat-np.array([[I(F(200,51))*z for z in row] for row in GB],dtype=object)
 for i in range(n):
  for j in range(i):KK[i,j]=KK[j,i]=I(max(KK[i,j].l,KK[j,i].l),min(KK[i,j].h,KK[j,i].h),True)
 test=certify(KK.tolist())
 print('parity',k,'exact projected band residual Gram row upper',float(upper(GB)),'corrected checked-band positive pivots',test['positive_pivots'],flush=True)
 kd=np.array([[float(mpiv(x)) for x in row] for row in KK]);evals,vecs=np.linalg.eigh(kd)
 print('parity',k,'checked-band correction min eigenvalue diagnostic',evals[0],flush=True)
 witness=None
 if evals[0]<0:
  z=[I(F(round(float(x)*10**12),10**12)) for x in vecs[:,0]]
  val=sum(z[i]*KK[i,j]*z[j] for i in range(n) for j in range(n));
  if val.h<0:witness=dict(z=[compact(x) for x in z],value=compact(val));print('exact negative projected finite comparison witness',float(F(val.h,S)),flush=True)
 out=dict(parity=k,cutoff=85,trial_band_end=KTRIAL,checked_end=N,finite_margin=str(m),finite_pivots=cert,checked_band_comparison=test,checked_band_gram_row_upper=str(upper(GB)),witness=witness,checked_band_matrix=[[compact(z) for z in row] for row in KK],scale=str(S))
 (ROOT/f'projected_{k}.json').write_text(json.dumps(out,indent=2)+'\n')
