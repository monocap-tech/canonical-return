"""Higher-order signed moment expansion for the energy-normalized tail.
Exact interval coefficients, no new moment constraints beyond BM's first four."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,sys
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2
root=Path(__file__).resolve().parent;p=int(sys.argv[1]);par=['even','odd'][p];W=4096;r=16;S=10**75
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 t=sum(int(v)**2 for v in a.flat);v=isqrt(t);return F(v+(v*v<t),s)
a=np.load(root/(par+'_corrector.npz'));j=a;C=obj(a['coefficient']);CS=int(a['coefficient_scale']);ec=F(str(a['coefficient_error']));nc=norm(C,CS);ns=[p]+list(map(int,j['indices']));w=[psi(n)[0] for n in ns];assert all(v.h<W*D for v in w);rr=[(F(1,4)+w[0]**2)/(F(1,4)+x*x) for x in w[1:]]
# Full coefficient l1 bound including pole reference.
A=(isqrt(len(ns)-1)+1+F(isqrt(sum(v.absmax()**2 for v in rr))+1,D))*(nc+ec)
pp={}
for q in range(2,55):
 if all(q%i for i in range(2,isqrt(q)+1)):
  k=q
  while k<=54:pp[k]=q;k*=q
assert len(pp)==24
wt=[(logi(I(q))/sqrtint(k),logi(I(k))) for k,q in pp.items()];d=[x*psi(n)[2]/2+sum((a*x*trig(x*s,1) for a,s in wt),I(0)) for n,x in zip(ns,w)]
aks=[F(4*k+1,2) for k in range(1,40)];E=[e2**(4*k+1) for k in range(1,40)];rat=[[(a*a)/(a*a+x*x) for x in w] for a in aks];a40=F(161,2)
# For every l<=r, ratio of successive exponential moments from k=40 <=1/2.
assert (e2**4*((a40+2)/a40)**(2*r)).h<D//2
# normalized absolute exponential moment sum <=1, used for remainder.
assert (sum((v*(a/W)**(2*r) for a,v in zip(aks,E)),I(0))+2*e2**161*(a40/W)**(2*r)).h<D

def rowbound(v):
 v=[x-t*v[0] for x,t in zip(v[1:],rr)];l=np.array([x.l*S//D for x in v],object);h=np.array([-((-x.h*S)//D) for x in v],object);m=(l+h)//2;er=norm(np.maximum(m-l,h-m),S);return norm(m@C,S*CS)+norm(m,S)*ec+er*(nc+ec)
ms=[];hs=[]
for l in range(2,r):
 m=[(x/W)**(2*l) for x in w];h=[a*b/W for a,b in zip(d,m)]
 for ak,ev,rv in zip(aks,E,rat):
  t=(-1)**(p+l)*ev*(ak/W)**(2*l)/W;h=[x+t*y for x,y in zip(h,rv)]
 err=(2*e2**161*(a40/W)**(2*l)/W).h;h=[v.widen(err) for v in h]
 ms.append(rowbound(m));hs.append(rowbound(h))
 if l in [2,3,4,8,16,24,32,39]:print(par,'l',l,'M normalized',float(ms[-1]),'H normalized',float(hs[-1]),flush=True)
# Since |d(x)|<=18x, nonexp remainder <=48 W^(2r) A/x^(2r+1) +48 W^(2r+1) A/x^(2r+2),
# exp remainder <=2 A W^(2r)/x^(2r+2). For x>=2W and W>=1, total <=74 W^(2r) A/x^(2r+1).
# Every term at cutoff decays at least as (N/m)^5. Integrate squared tail, with physical mass 2.
out={}
for N in [2**14,2**16,2**18,2**20,2**24]:
 x=F(3*N,4);assert x>=2*W
 head=sum((36*m*W**(2*l)/x**(2*l+1)+2*h*W**(2*l+1)/x**(2*l+2) for l,m,h in zip(range(2,r),ms,hs)),F(0));rem=74*W**(2*r)*A/x**(2*r+1);tau=F(N,18)*(head+rem)**2;crude=F(N,18)*(74*W**4*A/x**5)**2;assert tau<crude
 print(par,'N',N,'structured squared tail',float(tau),'crude',float(crude),'point remainder',float(rem),flush=True);out[str(N)]=str(tau)
np.savez_compressed(root/(par+'_U_tail.npz'),M=np.array(list(map(str,ms))),H=np.array(list(map(str,hs))),A=str(A),cutoff=np.array(list(map(int,out))),tail=np.array(list(out.values())),order=r)
print(par,'BO energy-normalized structured tail PASS',flush=True)
