"""Directed signed moment tail on the actual frozen four-column combined lift.
Uses inherited BM/BO moment cancellations and BO order-40 remainder, not
a claim of new moment-null constraints. Coefficient uncertainties retained.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,lcm
import numpy as np,sys,json
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'cr_int_t1d_ca'))
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2
S=10**75; W=4096; r=40; N=16384; G=10**12
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 t=sum(int(v)**2 for v in np.asarray(a).flat);v=isqrt(t);return F(v+(v*v<t),s)
def up(x,s=10**60):return F(-((-x.numerator*s)//x.denominator),s)
ca=root/'cr_int_t1d_ca';bo=root/'cr_int_t1d_bo'
p=np.load(ca/'bu_trial.npz');Q=obj(p['Q']);X=obj(p['X']);x=X@Q
j=np.load(bo/'even_joint.npz');b=np.load(bo/'even_normalization.npz');c=np.load(ca/'even_merged_corrector.npz')
B=obj(b['coefficient']);C=obj(c['coefficient']);bs=int(b['coefficient_scale']);cs=int(c['coefficient_scale']);eB=F(str(b['coefficient_error']));eC=F(str(c['coefficient_error']))
assert int(p['scale'])==G and Q.shape[1]==4 and B.shape==(len(j['indices']),Q.shape[0]) and C.shape==(len(c['indices']),x.shape[0])
assert eB>=0 and eC>=0
j1=j['indices'];j2=c['indices'];ids=np.union1d(j1,j2)
assert len(np.unique(j1))==len(j1) and len(np.unique(j2))==len(j2) and all(int(n)>0 and int(n)%2==0 for n in ids)
CS=lcm(bs*G,cs*G*G);Y=np.zeros((len(ids),4),object);U=Y.copy()
Y[np.searchsorted(ids,j1)]=(B@Q)*(CS//(bs*G))
U[np.searchsorted(ids,j2)]=(C@x)*(CS//(cs*G*G))
ey=[eB*norm(Q[:,k],G) for k in range(4)];eu=[eC*norm(x[:,k],G*G) for k in range(4)]
# The exact lift and exact correctors are separately pole-null and satisfy
# M0,H0,M1,H1=0 by inherited BO and BQ/CA constraints, hence so does their sum.
ns=[0]+list(map(int,ids));w=[psi(n)[0] for n in ns];assert all(v.h<W*D for v in w)
rr=[(F(1,4)+w[0]**2)/(F(1,4)+v*v) for v in w[1:]]
a_factor=isqrt(len(ns)-1)+1+F(isqrt(sum(v.absmax()**2 for v in rr))+1,D)
families={'Y':(Y,ey),'U':(U,eu),'combined':(Y+U,[a+b for a,b in zip(ey,eu)])}
norms={tag:[norm(M[:,k],CS)+ec[k] for k in range(4)] for tag,(M,ec) in families.items()}
A={tag:[a_factor*v for v in nn] for tag,nn in norms.items()}
pp={}
for q in range(2,55):
 if all(q%i for i in range(2,isqrt(q)+1)):
  k=q
  while k<=54:pp[k]=q;k*=q
assert len(pp)==24
wt=[(logi(I(q))/sqrtint(k),logi(I(k))) for k,q in pp.items()]
d=[v*psi(n)[2]/2+sum((a*v*trig(v*s,1) for a,s in wt),I(0)) for n,v in zip(ns,w)]
aks=[F(4*k+1,2) for k in range(1,40)];E=[e2**(4*k+1) for k in range(1,40)];rat=[[(a*a)/(a*a+v*v) for v in w] for a in aks];a40=F(161,2)
assert (e2**4*((a40+2)/a40)**(2*r)).h<D//2
assert (sum((v*(a/W)**(2*r) for a,v in zip(aks,E)),I(0))+2*e2**161*(a40/W)**(2*r)).h<D
def rowbound(v):
 v=[a-t*v[0] for a,t in zip(v[1:],rr)]
 l=np.array([a.l*S//D for a in v],object);h=np.array([-((-a.h*S)//D) for a in v],object);m=(l+h)//2
 er=norm(np.maximum(m-l,h-m),S);nm=norm(m,S)
 return {tag:[up(norm(np.array([val],object),S*CS)+nm*ec[k]+er*norms[tag][k]) for k,val in enumerate(m@M)] for tag,(M,ec) in families.items()}
ms=[];hs=[]
for l in range(2,r):
 m=[(v/W)**(2*l) for v in w];h=[a*b/W for a,b in zip(d,m)]
 for ak,ev,rv in zip(aks,E,rat):
  t=(-1)**l*ev*(ak/W)**(2*l)/W;h=[v+t*y for v,y in zip(h,rv)]
 err=(2*e2**161*(a40/W)**(2*l)/W).h;h=[v.widen(err) for v in h]
 ms.append(rowbound(m));hs.append(rowbound(h))
 print('CJ moment',l,{tag:[float(v) for v in ms[-1][tag]] for tag in families},flush=True)
xx=F(3*N,4);assert xx>=2*W
def sqru(x):
 s=10**45;k=-((-x.numerator*s*s)//x.denominator);t=isqrt(k);return F(t+(t*t<k),s)
result={}
for tag in families:
 tails=[]
 for k in range(4):
  h=sum((36*m[tag][k]*W**(2*l)/xx**(2*l+1)+2*v[tag][k]*W**(2*l+1)/xx**(2*l+2) for l,m,v in zip(range(2,r),ms,hs)),F(0))
  rem=74*W**(2*r)*A[tag][k]/xx**(2*r+1)
  tails.append(up(F(N,18)*(h+rem)**2))
 result[tag]={'squared':list(map(str,tails)),'norms':[str(sqru(v)) for v in tails],'frobenius':str(sqru(sum(tails)))}
 print('CJ TAIL',tag,[float(v) for v in tails],float(sum(tails)),flush=True)
out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True)
(out/'contracted-tail.json').write_text(json.dumps({'cutoff':N,'order':r,'families':result,'moments_M':[{k:list(map(str,v)) for k,v in q.items()} for q in ms],'moments_H':[{k:list(map(str,v)) for k,v in q.items()} for q in hs],'A':{k:list(map(str,v)) for k,v in A.items()},'coefficient_errors':{k:list(map(str,v[1])) for k,v in families.items()}},indent=2)+'\n')
