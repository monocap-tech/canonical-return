#!/usr/bin/env python3
"""CR-ALP-046: complete original signed finite E85, direct Gamma evaluation.
Classical coefficient and Euler-Maclaurin identities are analytical inputs.
No floating point in acceptance; no mixed Schur or Lean proof.
"""
from fractions import Fraction as F
from math import factorial,isqrt,comb
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
S=10**340; a=F(11,10)
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
M=90;SHIFT=512;N=85
assert F(4*factorial(2*M),6**(2*M)*2*M*SHIFT**(2*M))<F(1,10**290)
assert F(4*factorial(2*M),6**(2*M)*SHIFT**(2*M+1))<F(1,10**290)
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
 return re.widen(10**50),im.widen(10**50),tr.widen(10**50)
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
nu={j:j*NUF for j in range(1,N+1)};trig={j:[sincos(v*t) for t in TS] for j,v in nu.items()}
R={};R2={};G={};AI={}
for j,v in nu.items():
 re,im,tre=psi_pair(qt,v/2);R[j]=im/(2*v);R2[j]=im/(4*v*v*v)-tre/(8*v.sq())
 G[j]=re-log(PI);theta=sum(w*t[0] for w,t in zip(W,trig[j]));sgn=(-1)**j
 r=R[j]-sgn*sum(e/(b+v.sq()) for e,b in zip(EX,B2)).widen(1)
 R2[j]-=sgn*sum(e/(b+v.sq()).sq() for e,b in zip(EX,B2)).widen(1)
 AI[j]=v*r+theta
 if j%20==0:print('digamma complete',j,flush=True)
def Hss(i,j):
 u,v=nu[i],nu[j]
 if i==j:return G[i]-2*u.sq()/A*R2[i]-sum(w*((L-t)/A*c+s/(A*u)) for w,t,(s,c) in zip(W,TS,trig[i]))
 return -2/A*(v*AI[i]-u*AI[j])/(v.sq()-u.sq())
def integ(z):
 if z==0:return L
 return (E.sq()-1) if z==1 else (1-q.sq())
def integ_t(z):
 if z==0:return L.sq()/2
 assert z==-1
 return 1-(L+1)*q.sq()
def build(k):
 ids=list(range(k,N+1,2));sgn=(-1)**k;s=-sgn
 b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4
 Bmass=(E-q)/2+(A if k==1 else -A);ep=(b0+2*b1)/2;en=(b0-2*b1)/2
 E0=ep*L+en*(1-q.sq())-2*b0
 EC=[s*e*(b*b0-b1)/(bb-qt) for b,bb,e in zip(B[1:],B2[1:],EX[1:])]
 coef={F(1,2):ep,F(-1,2):en};kk=ZERO
 for p,cp in coef.items():
  for r,cr in coef.items():
   val=integ_t(p-F(1,2)) if r+F(1,2)==0 else integ(p+r)-integ(p-F(1,2))
   kk+=2*cp*cr*val
 j0=4*Bmass-kk
 jtail=-Bmass/2*S2+2*b0.sq()*S1+2*(b0.sq()/4-b1.sq())*S4-2*s*sum(e*(b*b0-b1).sq()/(bb-qt).sq() for b,bb,e in zip(B[1:],B2[1:],EX[1:])).widen(1)
 Pbb=ZERO
 for n,t,w in zip(NS,TS,W):
  sn=I(n).sqrt();ept={F(1,2):sn,F(-1,2):1/sn};z=ZERO
  for p,cp in coef.items():
   for r,cr in coef.items():
    v=p+r;it=L-t if v==0 else (E.sq()/n-1) if v==1 else (1-n*q.sq())
    z+=2*cp*cr*ept[r]*it
  Pbb-=w*z
 Hbb=g0*Bmass+j0+jtail+Pbb
 c={j:2*b0*nu[j]/(rt*(nu[j].sq()+qt)) for j in ids};Hb={}
 for j in ids:
  v=nu[j];den=v.sq()+qt;Lam=G[j]-g0;R1=R[j]-1/den;S3=Lam/(2*v.sq())-2/den
  Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v.sq()) for e,bb in zip(EC,B2[1:])).widen(1)
  val=g0*c[j]+Lam*(c[j]-b0/(rt*v))-2*v/rt*Es
  for n,w,(st,ct) in zip(NS,W,trig[j]):
   sn=I(n).sqrt();bt=ep*sn+en/sn;val-=2*w/(rt*den)*(v*(bt+b0*ct)-b1*st)
  Hb[j]=val
 # exact projected basis Pi e_j then congruence to w_j and u_k.
 H=[[Hss(i,j)-c[i]/Bmass*Hb[j]-c[j]/Bmass*Hb[i]+c[i]*c[j]/Bmass.sq()*Hbb for j in ids] for i in ids]
 Mass=[[I(int(i==j))-c[i]*c[j]/Bmass for j in ids] for i in ids]
 order=ids[1:]+[k];pos={j:i for i,j in enumerate(ids)};lam={j:c[j]/c[k] for j in ids[1:]}
 def convert(X,i,j):
  z=X[pos[i]][pos[j]]
  if i!=k:z-=lam[i]*X[0][pos[j]]
  if j!=k:z-=lam[j]*X[pos[i]][0]
  if i!=k and j!=k:z+=lam[i]*lam[j]*X[0][0]
  return z
 return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]
def bareiss_prefix(A):
 n=len(A);A=[row[:] for row in A];prev=1;out=[]
 for k in range(n):
  pivot=A[k][k];out.append(pivot)
  if pivot<=0:break
  for i in range(k+1,n):
   for j in range(i,n):
    z=A[i][j]*pivot-A[i][k]*A[k][j];r,rem=divmod(z,prev);assert rem==0;A[i][j]=A[j][i]=r
  prev=pivot
 return out
def cert(H,Mass,ix,sigma):
 A=[[H[i][j]-sigma*Mass[i][j] for j in ix] for i in ix];n=len(A)
 for i in range(n):
  for j in range(i):A[i][j]=A[j][i]=I(max(A[i][j].l,A[j][i].l),min(A[i][j].h,A[j][i].h),True)
 C=[[(A[i][j].l+A[i][j].h)//2 for j in range(n)] for i in range(n)]
 rad=max(sum(max(C[i][j]-A[i][j].l,A[i][j].h-C[i][j]) for j in range(n)) for i in range(n))
 for i in range(n):C[i][i]-=rad
 return bareiss_prefix(C),rad
if __name__=='__main__':
 bank=[];results=[]
 for k in (1,2):
  ids,H,Mass=build(k);print('built parity',k,'dimension',len(ids),flush=True)
  for X in (H,Mass):
   for i in range(len(ids)):
    for j in range(i):X[i][j]=X[j][i]=I(max(X[i][j].l,X[j][i].l),min(X[i][j].h,X[j][i].h),True)
  bank.append(dict(parity=k,ids=ids,H=[[[v.l,v.h] for v in r] for r in H],M=[[[v.l,v.h] for v in r] for r in Mass]))
  for cutoff in (14,24,40,66,85):
   ix=[i for i,j in enumerate(ids) if j<=cutoff];p,rad=cert(H,Mass,ix,F(1,10**110));count=sum(x>0 for x in p);print('cutoff',cutoff,'parity',k,'positive pivots',count,'/',len(ix),'radius units',rad,flush=True);results.append(dict(cutoff=cutoff,parity=k,positive_pivots=count,dimension=len(ix),sigma='1e-110',radius_units=rad))
  (ROOT/'CR-ALP-046_full_E85_intervals.json').write_text(json.dumps(dict(scale=S,blocks=bank))+'\n')
  (ROOT/'CR-ALP-046_result.json').write_text(json.dumps(results,indent=2)+'\n')
