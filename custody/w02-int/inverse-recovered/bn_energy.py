"""Directed integer intervals for the complete original Neumann Gamma-minus-prime matrix.
Analytic imports: AF Neumann identity, Euler-Maclaurin periodic Bernoulli remainder.
No floating point used in interval construction or sign acceptance.
"""
from fractions import Fraction as F
from math import isqrt,comb,factorial
from functools import lru_cache
from pathlib import Path
import numpy as np, sys,time
D=10**110
class I:
 __slots__=('l','h')
 def __init__(self,x=0,y=None,raw=False):
  if raw:self.l,self.h=int(x),int(y)
  else:
   a,b=F(x),F(x if y is None else y); self.l=a.numerator*D//a.denominator;self.h=-((-b.numerator*D)//b.denominator)
  assert self.l<=self.h
 def __add__(self,o):
  o=iv(o);return I(self.l+o.l,self.h+o.h,True)
 __radd__=__add__
 def __neg__(self):return I(-self.h,-self.l,True)
 def __sub__(self,o):return self+-iv(o)
 def __rsub__(self,o):return iv(o)+-self
 def __mul__(self,o):
  o=iv(o);v=[self.l*o.l,self.l*o.h,self.h*o.l,self.h*o.h];return I(min(v)//D,-((-max(v))//D),True)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=iv(o);assert o.l*o.h>0
  v=[(a*D,b) if b>0 else (-a*D,-b) for a in (self.l,self.h) for b in (o.l,o.h)]
  return I(min(a//b for a,b in v),max(-((-a)//b) for a,b in v),True)
 def __rtruediv__(self,o):return iv(o)/self
 def __pow__(self,n):
  assert n>=0
  r=I(1); a=self
  while n:
   if n&1:r=r*a
   a=a*a;n//=2
  return r
 def absmax(self):return max(abs(self.l),abs(self.h))
 def widen(self,e):return I(self.l-int(e),self.h+int(e),True)
 def __repr__(self):return str((self.l/D,self.h/D))
def iv(o):return o if isinstance(o,I) else I(o)
def sqrtint(n):
 r=isqrt(n*D*D);return I(r,r+1,True)
def atsmall(x):
 # |x|<=1/2, alternating power series, bound remainder by next absolute term.
 assert x.absmax()<=D//2
 t=x;s=t;x2=x*x
 for j in range(1,190):t=-t*x2;s+=t/(2*j+1)
 return s.widen((t*x2/381).absmax())
pi=16*atsmall(I(F(1,5)))-4*atsmall(I(F(1,239)))
def atanpos(x):
 assert x.l>=0
 if x.h<=D//2:return atsmall(x)
 if x.l>=2*D:return pi/2-atsmall(1/x)
 # atan(x)=pi/4+atan((x-1)/(x+1)); valid on [1/2,2].
 y=(x-1)/(x+1);assert y.absmax()<=D//2;return pi/4+atsmall(y)
def lognear(x):
 z=(x-1)/(x+1);assert z.absmax()<=D*34//100
 z2=z*z;t=z;s=t
 for k in range(1,140):t=t*z2;s+=t/(2*k+1)
 tail=(2*(t*z2)/281/(1-z2)).absmax();return (2*s).widen(tail)
log2=lognear(I(2))
def logi(x):
 assert x.l>0;k=0
 while x.l>=2*D:x=x/2;k+=1
 while x.h<D:x=x*2;k-=1
 return lognear(x)+k*log2
# e^-2 from alternating Taylor of exp(-2), terms decreasing after k=2.
t=I(1);ee=t
for k in range(1,160):t=-2*t/k;ee+=t
e2=ee.widen((2*t/160).absmax());assert e2.l>0 and e2.h<D//2
# complex intervals
class C:
 def __init__(self,a=0,b=0):self.a,self.b=iv(a),iv(b)
 def __add__(self,o):o=cv(o);return C(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return C(-self.a,-self.b)
 def __sub__(self,o):return self+-cv(o)
 def __mul__(self,o):o=cv(o);return C(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def inv(self):
  a,b=self.a,self.b;den=a*a+b*b;return C(a/den,-b/den)
 def __truediv__(self,o):return self*cv(o).inv()
def cv(o):return o if isinstance(o,C) else C(o)
bern=[F(1)]
for k in range(1,65):bern.append(-sum((F(comb(k+1,j))*bern[j] for j in range(k)),F(0))/F(k+1))
K=96;P=32;t0=F(K)+F(1,4)
# For Re z=t0, EM bound by integral of |B_{2P}(periodic)| |z+t|^{-2P-1}.
rempsi=I(abs(bern[2*P])/(2*P*t0**(2*P))).h
remtri=I(abs(bern[2*P])/t0**(2*P+1)).h
print('EM remainder bounds',rempsi,'/',D,remtri,'/',D,flush=True)
@lru_cache(None)
def psi(n):
 w=pi*int(n)/4;x=w/2;z=C(t0,x);q=z.inv();q2=q*q;power=q2
 a=C(logi(z.a*z.a+z.b*z.b)/2,atanpos(z.b/z.a))-q/2
 b=q+q2/2
 for k in range(1,P+1):
  a=a-F(bern[2*k],2*k)*power
  b=b+bern[2*k]*power*q;power=power*q2
 for k in range(K):
  r=C(F(k)+F(1,4),x).inv();a=a-r;b=b+r*r
 return w,a.a.widen(rempsi),a.b.widen(rempsi),b.a.widen(remtri)
# Complete nonexponential boundary sum as divided difference of w^2 S(w).
# S(w)=sum 1/(a_k^2+w^2)=Im psi(.25+iw/2)/(2w) for w>0.
# At zero S=psi1(.25)/4. Bdiag=2S+wS'=S+Re psi1/4 for w>0.
def trig(x,odd):
 # Reduce using rational midpoint to nearest 2pi; interval carries uncertainty.
 k=(x.l+x.h+(pi.l+pi.h))//(2*(pi.l+pi.h));x=x-2*k*pi;assert x.absmax()<4*D
 t=x if odd else I(1); s=t;x2=x*x
 for k in range(1,90):t=-t*x2/((2*k+odd-1)*(2*k+odd));s+=t
 return s.widen((t*x2/((180+odd-1)*(180+odd))).absmax())
def build(ns,par):
 w=[];m=[];S=[];Bdiag=[];J=[]
 for n in ns:
  x,a,b,c=psi(int(n));s=b/(2*x) if n else c/4
  w.append(x);m.append(a-logi(pi));S.append(s);Bdiag.append(s+c/4 if n else 2*s);J.append(x*x*s)
 w2=[x*x for x in w];N=len(ns);A=[[I(0) for _ in ns] for _ in ns];print('digamma scalar intervals done',N,flush=True)
 pp={}
 for p in range(2,55):
  if all(p%k for k in range(2,isqrt(p)+1)):
   n=p
   while n<=54:pp[n]=p;n*=p
 assert len(pp)==24
 for n,p in pp.items():
  s=logi(I(n));wt=logi(I(p))/sqrtint(n);sn=[trig(x*s,1) for x in w];cs=[trig(x*s,0) for x in w]
  for i in range(N):
   for j in range(i+1):
    if i==j:v=2*(4-s) if ns[i]==0 else (4-s)*cs[i]-sn[i]/w[i]
    else:v=(sn[j]-sn[i])/(w[i]-w[j])-(sn[i]+sn[j])/(w[i]+w[j])
    A[i][j]+=wt*v
 print('all 24 signed arithmetic matrices done',flush=True)
 # exponential correction after pole-null congruence. Use rank-one terms with k>=1;
 # k=0 vanishes exactly. Nonexp full B includes k=0, annihilated on v_n.
 r=[(F(1,4)+w2[0])/(F(1,4)+x) for x in w2[1:]]
 mass=[4 if n==0 else 2 for n in ns]
 G=[[I(0) for _ in ns] for _ in ns]
 for i in range(N):
  for j in range(i+1):
   v=Bdiag[i] if i==j else 2*(J[i]-J[j])/(w2[i]-w2[j]);G[i][j]=G[j][i]=v-A[i][j]+(mass[i]*m[i] if i==j else 0)
 out=[[G[i+1][j+1]-r[i]*G[0][j+1]-r[j]*G[i+1][0]+r[i]*r[j]*G[0][0] for j in range(N-1)] for i in range(N-1)]
 # exp tail k>=40: since a_k>=2k and exp(-4a_k)<=e^(-2)^(4k+1),
 # each projected boundary vector <=(1+|r|)/a_k^2<=2/a_k^2.
 # absolute entry <=8 e^(-4a_k)/a_k^2 <=8 e^(-4a_k), geometric sum.
 for k in range(1,40):
  a=F(2*k)+F(1,2);t=[1/(a*a+w2[i+1])-r[i]/(a*a+w2[0]) for i in range(N-1)]
  c=(-1)**(par+1)*2*a*a*(e2**(4*k+1))
  for i in range(N-1):
   for j in range(i+1):out[i][j]+=c*t[i]*t[j]
 assert all(x.l>=0 and x.h<=D for x in r)
 tail=(8*e2**161/(1-e2**4)).h
 for i in range(N-1):
  for j in range(i+1):out[i][j]=out[i][j].widen(tail);out[j][i]=out[i][j]
 return out,r
if __name__=='__main__':
 p=int(sys.argv[1]);par=['even','odd'][p];root=Path(__file__).resolve().parent;suffix='_extension' if '--extend' in sys.argv else '';ar=np.load(root/(par+suffix+'_lift.npz'));ns=[p]+list(map(int,ar['indices']));start=time.time();Q,r=build(ns,p)
 width=max(x.h-x.l for row in Q for x in row);print('complete signed interval width',width,'/',D,'diagnostic',width/D,'time',time.time()-start,flush=True)
 np.savez_compressed(root/(par+suffix+'_energy.npz'),lo=np.array([[str(x.l) for x in row] for row in Q]),hi=np.array([[str(x.h) for x in row] for row in Q]),scale=str(D),indices=ns[1:])
 # Direct interval LDL on full selected lift matrix (not Schur). No heuristic sign accepted.
 A=[row[:] for row in Q];piv=[]
 for i in range(len(A)):
  pivot=A[i][i];piv.append(pivot)
  if pivot.l<=0:break
  v=[A[j][i]/pivot for j in range(i+1,len(A))]
  for j in range(i+1,len(A)):
   for k in range(j,len(A)):A[k][j]=A[j][k]=A[k][j]-A[k][i]*v[j-i-1]
 print('interval LDL positive pivots',sum(x.l>0 for x in piv),'/',len(A),'last',piv[-1],flush=True)
 if len(piv)==len(A) and all(x.l>0 for x in piv):print('BL signed finite lift compression positive PASS',par,flush=True)
 else:print('BL lift sign NOT certified; no inference of negativity',par,flush=True)
 np.savez_compressed(root/(par+suffix+'_pivots.npz'),lo=np.array([str(x.l) for x in piv]),hi=np.array([str(x.h) for x in piv]),scale=str(D))
