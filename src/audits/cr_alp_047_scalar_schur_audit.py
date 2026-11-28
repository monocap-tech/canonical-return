#!/usr/bin/env python3
"""CR-ALP-047: finite-band scalar Schur falsifier; original signed coefficients.
Classical coefficient and Euler-Maclaurin identities are analytical inputs.
Midpoint search only; exact intervals verify selected witnesses. No Lean proof.
"""
from fractions import Fraction as F
from math import factorial,isqrt,comb
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
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
assert F(4)**800/factorial(800)<F(1,S)
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
 for m in range(1,400):
  ts=ts*q/((2*m)*(2*m+1));tc=tc*q/((2*m-1)*(2*m));s+=ts;c+=tc
 # cos through degree 798, sin through degree 799. Error <= 4^800/800!.
 return s.widen(1),c.widen(1)
# Euler-Maclaurin with derivative, shift 512 and 90 Bernoulli terms.
M=24;SHIFT=256;N=105
assert F(4*factorial(2*M),6**(2*M)*2*M*SHIFT**(2*M))<F(1,10**92)
assert F(4*factorial(2*M),6**(2*M)*SHIFT**(2*M+1))<F(1,10**92)
BER=[F(1)]
for n in range(1,2*M+1):BER.append(-sum(F(comb(n+1,k))*BER[k] for k in range(n))/(n+1))
def cmul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def inv(x,y):den=x.sq()+y.sq();return x/den,-y/den
def psi_pair(x,y=ZERO):
 x,y=iv(x),iv(y);xx=x+SHIFT;assert y.up()/F(xx.l,S)<F(1,3)
 lr=log(xx.sq()+y.sq())/2;li=atan(y/xx);z=inv(xx,y);u=cmul(z,z);p=u
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
import mpmath as mp
mp.mp.dps=85
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
def mpiv(x):return mp.mpf(x.l+x.h)/(2*S)
def fwd(C,r):
 z=mp.matrix(len(r),1)
 for i in range(len(r)):z[i]=(r[i]-sum(C[i,h]*z[h] for h in range(i)))/C[i,i]
 return z
def back(C,r):
 z=mp.matrix(len(r),1)
 for i in reversed(range(len(r))):z[i]=(r[i]-sum(C[h,i]*z[h] for h in range(i+1,len(r))))/C[i,i]
 return z
def rup(x):return str(F(ceildiv(x.h*10**50,S),10**50))
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
results=[]
for bl in bank['blocks']:
 k=bl['parity'];ids=bl['ids'];n=len(ids);sgn=(-1)**k;s=-sgn
 H=[[I(l//div,ceildiv(h,div),True) for l,h in row] for row in bl['H']]
 C=mp.cholesky(mp.matrix([[mpiv(x) for x in row] for row in H])); print('parity',k,'min energy pivot',mp.nstr(min(C[i,i]**2 for i in range(n)),10),flush=True)
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
 js=list(range(87 if k==1 else 86,N+1,2));rows={j:[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-c(k)/BM*hb(j) for i in ids] for j in js}
 Z={j:fwd(C,mp.matrix([mpiv(x) for x in rows[j]])) for j in js};best=None
 for j in js:
  for l in js:
   if l<=j:continue
   lam=mpiv(c(j)/c(l));z=Z[j]-lam*Z[l];rho=(z.T*z)[0]/(1+lam**2)
   if best is None or rho>best[0]:best=(rho,j,l,z)
 rho,j,l,z=best;p=back(C,z);p/=max(abs(x) for x in p);den=10**55;nums=[int(mp.nint(t*den)) for t in p];pp=[I(F(t,den)) for t in nums]
 # Exact verification, not relying on midpoint Cholesky or optimization.
 lam=c(j)/c(l); norm=1+lam.sq();cross=sum(x*(u-lam*v) for x,u,v in zip(pp,rows[j],rows[l]));energy=sum(pp[i]*H[i][h]*pp[h] for i in range(n) for h in range(n));assert energy.l>0
 ratio=cross.sq()/(norm*energy);print('pair',j,l,'ratio',float(F(ratio.l,S)),float(F(ratio.h,S)),'energy',float(F(energy.l,S)),flush=True)
 # certify strict failure of H-4 B when 4 cross² > W(p) ||q||²
 fail=4*cross.sq()-energy*norm;
 if k==1:assert ratio.l*100>53*S and fail.l>0
 qenergy=hss(j,j)-2*lam*hss(j,l)+lam.sq()*hss(l,l);assert qenergy.l>0
 det=energy*qenergy-cross.sq();print('actual q energy per norm',float((qenergy/norm).mid()),'actual 2d normalized coupling',float((cross.sq()/(energy*qenergy)).mid()),'actual two-dimensional positive',det.l>0,flush=True)
 print('scalar 1/4 sufficient gate falsified?',fail.l>0,flush=True)
 results.append(dict(parity=k,ids=ids,j=j,l=l,coefficient_numerators=list(map(str,nums)),coefficient_denominator=str(den),energy_interval=list(map(str,(energy.l,energy.h))),cross_interval=list(map(str,(cross.l,cross.h))),q_norm_sq_interval=list(map(str,(norm.l,norm.h))),ratio_interval=list(map(str,(ratio.l,ratio.h))),four_cross_sq_minus_energy_norm_interval=list(map(str,(fail.l,fail.h))),scalar_gate_falsified=fail.l>0,actual_q_energy_interval=list(map(str,(qenergy.l,qenergy.h))),actual_two_dim_det_interval=list(map(str,(det.l,det.h)))))
 # Multi-mode finite-band candidate, then exact pole-null lift.
 pole=mp.matrix([mpiv(c(jj)) for jj in js]); Zm=mp.matrix([[Z[jj][i] for i in range(n)] for jj in js]); PZ=Zm-pole*(pole.T*Zm)/(pole.T*pole)[0]
 eig,V=mp.eigsy(PZ*PZ.T); y=V[:,len(js)-1]; y/=max(abs(t) for t in y)
 yn=[int(mp.nint(t*den)) for t in y[:-1]]; yc=[I(F(t,den)) for t in yn];yc.append(-sum(c(jj)*t for jj,t in zip(js[:-1],yc))/c(js[-1]))
 zr=sum((Z[jj]*mpiv(t) for jj,t in zip(js,yc)),mp.matrix(n,1)); pv=back(C,zr);pv/=max(abs(t) for t in pv);pn=[int(mp.nint(t*den)) for t in pv];pc=[I(F(t,den)) for t in pn]
 pe=sum(pc[i]*H[i][h]*pc[h] for i in range(n) for h in range(n)); assert pe.l>0
 qn=sum(t.sq() for t in yc);assert qn.l>0
 cr=sum(t*sum(pc[i]*rows[jj][i] for i in range(n)) for jj,t in zip(js,yc));qe=sum(yc[i]*hss(jj,ll)*yc[h] for i,jj in enumerate(js) for h,ll in enumerate(js));assert qe.l>0
 rr=cr.sq()/(pe*qn); dd=pe*qe-cr.sq();actual=cr.sq()/(pe*qe)
 print('multi-band parity',k,'ratio',float(rr.mid()),'actual q energy/norm',float((qe/qn).mid()),'actual coupling fraction',float(actual.mid()),'actual determinant positive',dd.l>0,flush=True)
 def pair(x):return list(map(str,(x.l,x.h)))
 results[-1]['multi_band']=dict(js=js,p_numerators=list(map(str,pn)),q_first_numerators=list(map(str,yn)),denominator=str(den),ratio_interval=pair(rr),actual_coupling_fraction_interval=pair(actual),q_energy_per_norm_interval=pair(qe/qn),p_energy_interval=pair(pe),q_energy_interval=pair(qe),q_norm_sq_interval=pair(qn),cross_interval=pair(cr),actual_det_interval=pair(dd))
 # Entire finite pole-null band, not just individual witnesses.
 assert dd.l>0 and det.l>0
 if k==1:assert rr.l*100>105*S
 else:assert rr.l*100>27*S
 last=js[-1];js0=js[:-1];ls=[c(jj)/c(last) for jj in js0];nb=len(js0)
 Bcross=[[rows[jj][i]-la*rows[last][i] for jj,la in zip(js0,ls)] for i in range(n)]
 HH=[[hss(jj,ll)-lb*hss(jj,last)-la*hss(last,ll)+la*lb*hss(last,last) for ll,lb in zip(js0,ls)] for jj,la in zip(js0,ls)]
 # Symmetry by explicit intersection of independent interval expressions.
 full=[H[i]+Bcross[i] for i in range(n)]+[[Bcross[i][h] for i in range(n)]+HH[h] for h in range(nb)]
 for i in range(n+nb):
  for j in range(i):full[i][j]=full[j][i]=I(max(full[i][j].l,full[j][i].l),min(full[i][j].h,full[j][i].h),True)
 cert=certify(full);print('complete finite mixed band',k,cert,flush=True);assert cert['all_positive'];results[-1]['finite_band_joint_certificate']=cert
 (ROOT/'result.json').write_text(json.dumps(dict(scale=str(S),bank_sha256=hashlib.sha256(raw).hexdigest(),blocks=results),indent=2)+'\n')
