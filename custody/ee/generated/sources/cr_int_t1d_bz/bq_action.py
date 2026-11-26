"""Tiled exact-grid signed action in BO normalized coordinates. Inherits BO source identities.
Nonexp divided differences + separable exp correction; exact integer matrix products.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import sys,time,numpy as np
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2,C,bern
root=Path(__file__).resolve().parent;bo=root.parent/'cr_int_t1d_bo';p=int(sys.argv[1]);par=['even','odd'][p];cut=int(sys.argv[2]) if len(sys.argv)>2 else 16384
S=10**70;BS=10**65;YS=10**12
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 k=sum(int(v)**2 for v in a.flat);r=isqrt(k);return F(r+(r*r<k),s)
def grid(rows):
 l=np.array([[v.l*S//D for v in row] for row in rows],object);h=np.array([[-((-v.h*S)//D) for v in row] for row in rows],object);a=(l+h)//2;return a,norm(np.maximum(a-l,h-a),S)
def up(x,s=10**45):return F(-((-x.numerator*s)//x.denominator),s)
a=np.load(root/(par+'_corrector.npz'));ns=[p]+list(map(int,a['indices']));lo=int(sys.argv[3]) if len(sys.argv)>3 else p;ms=list(range(lo,cut+1,2));assert lo%2==p;w={n:pi*n/4 for n in set(ns+ms)};w2={n:v*v for n,v in w.items()};ref=w2[p];rr=[(F(1,4)+ref)/(F(1,4)+w2[n]) for n in ns[1:]];bb=obj(a['coefficient']);bs=int(a['coefficient_scale']);B=(bb*BS)//bs;eB=up(F(str(a['coefficient_error']))+F(len(B.flat),BS));nB=norm(B,BS)
# Direct EM at .25+i*x, no 96-term shift needed when m>256. Integral of |z+t|^(-2P-1)
# <= x^(-2P)*(1+1/(2P)), using |z+t|>=max(x,t). Periodic Bernoulli remainder <= |B_2P| times integral.
def fast_im(n):
 if n<=256 or n in ns:return psi(n)[2]
 x=w[n]/2;assert x.l>100*D;v=F(1,4)/x;t=v;ar=t
 for k in range(1,28):t=-t*v*v;ar+=t/(2*k+1)
 ar=ar.widen((t*v*v/57).absmax());z=C(F(1,4),x);q=z.inv();q2=q*q;pow=q2;res=pi/2-ar-q.b/2
 for k in range(1,33):res-=F(bern[2*k],2*k)*pow.b;pow=pow*q2
 rem=I(abs(bern[64])*F(65,64)/(F(x.l,D)**64)).h
 return res.widen(rem)
# all 24 prime power pairs
pp={}
for q in range(2,55):
 if all(q%i for i in range(2,isqrt(q)+1)):
  k=q
  while k<=54:pp[k]=q;k*=q
assert len(pp)==24
wt=[(logi(I(q))/sqrtint(k),logi(I(k))) for k,q in pp.items()]
# Direct phases for selected modes. For consecutive head rows, exact interval complex recurrence,
# reseeded every 32 steps to limit interval width. No floating point phase propagation.
sn={};cs={}
for n in ns:
 sn[n]=[trig(w[n]*t,1) for _,t in wt];cs[n]=[trig(w[n]*t,0) for _,t in wt]
d={n:w[n]*psi(n)[2]/2+w[n]*sum((a*s for (a,_),s in zip(wt,sn[n])),I(0)) for n in ns}
mG={n:psi(n)[1]-logi(pi) for n in ns};bd={n:(psi(n)[2]/(2*w[n])+psi(n)[3]/4 if n else psi(n)[3]/2) for n in ns};ad={n:sum((a*(2*(4-t) if n==0 else (4-t)*c-s/w[n]) for (a,t),s,c in zip(wt,sn[n],cs[n])),I(0)) for n in ns}
# Independent fast-EM/direct check (direct psi uses shifted EM).
for n in [258+p,510+p,1022+p,4094+p,8190+p,16382+p]:
 if lo<=n<=cut:
  v=fast_im(n);v0=psi(n)[2];assert max(v.l,v0.l)<=min(v.h,v0.h),(par,n)
print(par,'scalar prelude and EM overlap PASS',flush=True)
step=[(trig(pi*t/2,0),trig(pi*t/2,1)) for _,t in wt];mass=sqrtint(2);fac=32/(pi*pi*mass)
# Factor exponential sum k=1..39; reference pole k=0 vanishes exactly.
ex=[(F(4*k+1,2),(-1)**(p+1)*2*F(4*k+1,2)**2*e2**(4*k+1)) for k in range(1,40)];tr=[[1/(ak*ak+w2[n])-r/(ak*ak+ref) for n,r in zip(ns[1:],rr)] for ak,_ in ex];T,eT=grid(tr);nT=norm(T,S);TB=T@B;tail=F((4*e2**161/(1-e2**4)).h,D)
# For off-diagonal raw/mass: common fac/(m^2-n^2). Avoid interval square subtraction.
# For m>256 base reference never coincides.
Ytiles=[];tile_err=[];start=time.time();new=ms;tests=0
for begin in range(0,len(new),256):
 tile=new[begin:begin+256];rows=[];left=[]
 for i,m in enumerate(tile):
  if i%32==0:phase=[(trig(w[m]*t,0),trig(w[m]*t,1)) for _,t in wt]
  else:phase=[(a*c-b*s,a*s+b*c) for (a,b),(c,s) in zip(phase,step)]
  dm=d[m] if m in d else w[m]*fast_im(m)/2+w[m]*sum((a*s for (a,_),(_,s) in zip(wt,phase)),I(0))
  den=I(2) if m==0 else mass
  # all rows, including the pole reference and diagonals
  def entry(n):
   if m==n:return (bd[n]-ad[n]+(4 if m==0 else 2)*mG[n])/den
   return 2*(dm-d[n])/(w2[m]-w2[n])/den
  r0=entry(p);row=[entry(n)-r*r0 for n,r in zip(ns[1:],rr)]
  rows.append(row);left.append([coef/(ak*ak+w2[m])/den for ak,coef in ex])
  # independent direct phase checks on selected tile endpoints
  if i==len(tile)-1:
   for (_,t),(c,s) in zip(wt,phase):
    for z,zz in [(c,trig(w[m]*t,0)),(s,trig(w[m]*t,1))]:assert max(z.l,zz.l)<=min(z.h,zz.h);tests+=1
 R,eR=grid(rows);X,eX=grid(left);nR=norm(R,S);nX=norm(X,S)
 # Exact signed contraction; central action R+XT. Remainder enclosed entrywise as BO.
 Y=(R@B*S+X@TB)*YS//(S*S*BS)
 eta=up((nR+nX*nT)*eB+(eR+eX*nT+(nX+eX)*eT+tail*len(tile)*len(rr))*(nB+eB)+F(Y.size,YS))
 assert max(abs(v) for v in Y.flat)<2**63
 Ytiles.append(Y.astype(np.int64));tile_err.append(eta)
 if begin%1024==0 or begin+len(tile)==len(new):print(par,'through',tile[-1],'tile error',float(eta),'elapsed',round(time.time()-start,1),flush=True)
 # compact partial checkpoints permit a valid smaller head if execution interrupted
 np.savez_compressed(root/(par+f'_U_{lo}_{cut}_partial.npz'),Y=np.concatenate(Ytiles),indices=ms[:sum(len(v) for v in Ytiles)],scale=YS,tile_errors=np.array(list(map(str,tile_err))))
Y=np.concatenate(Ytiles);assert len(Y)==len(ms);err=sum(tile_err,F(0));assert err<F(1,10**5)
print(par,'complete normalized head',cut,'rows',len(Y),'error',float(err),'Frobenius',float(norm(Y,YS)+err),'phase overlap checks',tests,flush=True)
np.savez_compressed(root/(par+f'_U_{lo}_{cut}.npz'),Y=Y,indices=ms,scale=YS,error=str(err),tile_errors=np.array(list(map(str,tile_err))))
print(par,'BP full head PASS',flush=True)
