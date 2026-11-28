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
M=24;SHIFT=1024;N=850;KTRIAL=600
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
 def pb(j):
  v=nu[j];den=v.sq()+qt;val=ZERO
  for nn,w,(st,ct) in zip(NS,W,trig[j]):
   sn=I(nn).sqrt();bt=ep*sn+en/sn;val-=2*w/(rt*den)*(v*(bt+b0*ct)-b1*st)
  return val
 def hb(j):
  v=nu[j];den=v.sq()+qt;Lam=G[j]-g0;R1=R[j]-1/den;S3=Lam/(2*v.sq())-2/den
  Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v.sq()) for e,bb in zip(EC,B2[1:])).widen(1)
  val=g0*c(j)+Lam*(c(j)-b0/(rt*v))-2*v/rt*Es
  for nn,w,(st,ct) in zip(NS,W,trig[j]):
   sn=I(nn).sqrt();bt=ep*sn+en/sn;val-=2*w/(rt*den)*(v*(bt+b0*ct)-b1*st)
  return val
 def load(A):return mat([[I(int(l),int(h),True) for l,h in row] for row in A])
 def midfloat(A):return np.array([[float(x.mid()) for x in row] for row in A])
 def rr(A):return mat([[I(F(round(float(x)*10**16),10**16)) for x in row] for row in A])
 def sym(A):
  for i in range(len(A)):
   for j in range(i): A[i,j]=A[j,i]=I(max(A[i,j].l,A[j,i].l),min(A[i,j].h,A[j,i].h),True)
  return A
 def save(A):return [[compact(x) for x in row] for row in A]
 # Exact integer midpoint/radius matrix product; substantially cheaper than interval-object products.
 def mul(A,B):
  def mr(A):
   M=np.array([[(x.l+x.h)//2 for x in row] for row in A],dtype=object)
   R=np.array([[max(m-x.l,x.h-m) for x,m in zip(row,ms)] for row,ms in zip(A,M)],dtype=object)
   return M,R
  ma,ra=mr(A);mb,rb=mr(B); mc=ma@mb; rc=abs(ma)@rb+ra@abs(mb)+ra@rb
  return mat([[I((x-r)//S,ceildiv(x+r,S),True) for x,r in zip(row,rad)] for row,rad in zip(mc,rc)])
 oldraw=(ROOT.parent/'cr_alp_049'/f'result_{k}.json').read_bytes();assert hashlib.sha256(oldraw).hexdigest()=={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k];old=json.loads(oldraw);V=load(old['V']);EV=sym(mul(mul(V.T,H),V))
 assert upper(EV-mat([[I(int(i==j)) for j in range(n)] for i in range(n)]))<F(1,10**12)
 band=list(range(87 if k==1 else 86,N+1,2));trial=[j for j in band if j<=KTRIAL];nt=len(trial)
 rows=mat([[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-c(k)/BM*hb(j) for i in ids] for j in band]);U=mul(rows,V)
 Lraw=mat([[(I(gamma) if i==j else ZERO)+pss(i,j) for j in trial] for i in band]);Ltt=sym(Lraw[:nt,:].copy());C=mat([[c(j)] for j in band]);ct=C[:nt,:]
 print('parity',k,'coefficient matrices built',flush=True)
 # g = (c_A - r R85 b)/sqrt(s_A(1-r)), r=s_A/s_tot.
 low=list(range(k,86,2)); allids=low+trial; cl=mat([[c(j)] for j in low]); ca=ct
 sl=mul(cl.T,cl)[0,0]; sa=mul(ca.T,ca)[0,0]; st=BM-sl; si=st-sa; r=sa/st; mass=sa*si/st; norm=mass.sqrt()
 assert si.l>0 and mass.l>0
 d=mat([[r*c(j)/norm] for j in low]+[[c(j)/norm] for j in trial]); db=-r/norm
 # W(b,b) recovered from the pinned original anchor diagonal.
 anchor=ids.index(k); alpha=c(k)/BM
 hbb=(H[anchor,anchor]-hss(k,k)+2*alpha*hb(k))/alpha.sq()
 hlow=mat([[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-alpha*hb(j) for i in ids] for j in low])
 hbp=mat([[hb(i)-c(i)/c(k)*hb(k) if i!=k else hb(k)-alpha*hbb for i in ids]])
 Bg=mul(d.T,np.concatenate((mul(hlow,V),U[:nt,:]),axis=0))+mul(hbp,V)*db
 # Exact original prime pole-pole coefficient from CR-ALP-046.
 Pbb=ZERO;coef={F(1,2):ep,F(-1,2):en}
 for nn,t,w in zip(NS,TS,W):
  sn=I(nn).sqrt();ept={F(1,2):sn,F(-1,2):1/sn};zz=ZERO
  for pp,cp in coef.items():
   for rr0,cr in coef.items():
    v=pp+rr0;it=L-t if v==0 else (E.sq()/nn-1) if v==1 else (1-nn*q.sq())
    zz+=2*cp*cr*ept[rr0]*it
  Pbb-=w*zz
 Pall=mat([[pss(i,j) for j in allids] for i in allids]); pball=mat([[pb(j)] for j in allids]); pg=mul(Pall,d)+pball*db
 Lgg=I(gamma)+(mul(d.T,pg)[0,0]+db*(mul(pball.T,d)[0,0]+db*Pbb))
 Lfg=pg[len(low):,:] # gamma bridge is pole-aligned on the finite band, vanishes against F.
 Uaug=np.concatenate((U[:nt,:],Bg),axis=0)
 Laug=np.block([[Ltt,Lfg],[Lfg.T,mat([[Lgg]])]])
 cf=np.concatenate((ct,mat([[ZERO]])),axis=0)
 cfmid=midfloat(cf); lm=midfloat(Laug); um=midfloat(Uaug)
 sad=np.block([[lm,cfmid],[cfmid.T,np.zeros((1,1))]])
 sol=np.linalg.solve(sad,np.vstack((um,np.zeros((1,n)))))[:-1,:]
 X=rr(sol);X[nt-1,:]=-mul(mat([[c(j)/c(trial[-1]) for j in trial[:-1]]]),X[:nt-1,:])[0,:]
 cross=mul(Uaug.T,X);Knew=sym(EV-cross-cross.T+mul(X.T,mul(Laug,X)))
 Raug=Uaug-mul(Laug,X); fa=Raug[:nt,:]; fa=fa-mul(ct,mat([[x/sa for x in mul(ct.T,fa)[0]]]))
 Gtrial=sym(mul(fa.T,fa)+mul(Raug[nt:,:].T,Raug[nt:,:]));
 print('parity',k,'bridge norm2 unnormalized',float(mass.mid()),'Lgg',float(Lgg.mid()),'trial residual upper',float(upper(Gtrial)),flush=True)
 Kfinite=sym(Knew-Gtrial*I(F(200,51)))
 for sh in [F(1,4),F(1,8),F(1,16),F(1,32),F(1,100),F(0)]:
  test=certify(Kfinite-mat([[I(sh if i==j else 0) for j in range(n)] for i in range(n)]))
  if test['all_positive']:break
 print('parity',k,'augmented finite test',str(sh),test,flush=True)
 high=band[nt:]; ch=C[nt:,:]
 Pgout=mul(mat([[pss(i,j) for j in allids] for i in high]),d)+mat([[pb(j)] for j in high])*db
 Lgout=Pgout+ch*(I(gamma)*db)
 Y=U[nt:,:]-mul(Lraw[nt:,:],X[:nt,:])-mul(Lgout,X[nt:,:])
 cy=mul(ch.T,Y); cc=mul(ch.T,ch)[0,0]; Yp=Y-mul(ch,mat([[x/cc for x in cy[0]]]))
 Gout=sym(mul(Yp.T,Yp));Cmp=sym(Kfinite-Gout*I(F(200,51)))
 ge,gv=np.linalg.eigh(midfloat(Gout));ce,cv=np.linalg.eigh(midfloat(Cmp))
 def witness(A,v):
  zz=rr(v.reshape(-1,1)); f=mul(zz.T,mul(A,zz))[0,0]/mul(zz.T,zz)[0,0];return {'z':save(zz),'lo':str(F(f.l,S)),'hi':str(F(f.h,S))}
 w=witness(Gout,gv[:,-1]);wc=witness(Cmp,cv[:,0]);print('parity',k,'outside projected norm2 witness',float(F(w['lo'])),'scalar comparison witness upper',float(F(wc['hi'])),flush=True)
 out=dict(parity=k,scale=str(S),trial=600,band_end=N,bridge_mass=compact(mass),bridge_energy=compact(Lgg),Knew=save(Knew),Gtrial=save(Gtrial),Kfinite=save(Kfinite),Gout=save(Gout),comparison=save(Cmp),finite_shift=str(sh),finite_test=test,norm_witness=w,comparison_witness=wc,trial_free_rows=save(X[:-2,:]),bridge_coeff=save(X[-1:,:]))
 (ROOT/f'result_{k}.json').write_text(json.dumps(out)+'\n')
