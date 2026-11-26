"""Rectangular original signed-action pilot; exact Frobenius projection upper bound for finite head.
Not an energy-normalized or whole-tail residual test."""
from pathlib import Path
import sys,numpy as np
from fractions import Fraction as F
from math import isqrt,factorial
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2
root=Path(__file__).resolve().parent;bk=root.parent/'cr_int_t1d_bk';p=int(sys.argv[1]);par=['even','odd'][p];cut=int(sys.argv[2]) if len(sys.argv)>2 else 256
S=10**75; obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 k=sum(int(v)**2 for v in a.flat);r=isqrt(k);return F(r+(r*r<k),s)
def rnd(a):return np.vectorize(lambda v:int(F(float(v))*S),otypes=[object])(a)
def fl(a,s=S):return np.vectorize(lambda v:float(F(int(v),s)),otypes=[float])(a)
def grid(rows):
 l=np.array([[v.l*S//D for v in row] for row in rows],object);h=np.array([[-((-v.h*S)//D) for v in row] for row in rows],object);a=(l+h)//2;return a,norm(np.maximum(a-l,h-a),S)
a=np.load(root/(par+'_joint.npz'));ns=[p]+list(map(int,a['indices']));ms=list(range(p,cut+1,2));alln=sorted(set(ns+ms));w={n:psi(n)[0] for n in alln};J={n:w[n]*psi(n)[2]/2 for n in alln};bd={n:(psi(n)[2]/(2*w[n])+psi(n)[3]/4 if n else psi(n)[3]/2) for n in alln};mG={n:psi(n)[1]-logi(pi) for n in alln};pp={}
for q in range(2,55):
 if all(q%i for i in range(2,isqrt(q)+1)):
  k=q
  while k<=54:pp[k]=q;k*=q
assert len(pp)==24
wt=[(logi(I(q))/sqrtint(k),logi(I(k))) for k,q in pp.items()];sn=[{n:trig(w[n]*t,1) for n in alln} for _,t in wt];cs=[{n:trig(w[n]*t,0) for n in alln} for _,t in wt];ad={n:sum((a*(2*(4-t) if n==0 else (4-t)*c[n]-s[n]/w[n]) for (a,t),s,c in zip(wt,sn,cs)),I(0)) for n in alln};d={n:J[n]+w[n]*sum((a*s[n] for (a,_),s in zip(wt,sn)),I(0)) for n in alln}
# Same signed off-diagonal identity as BM, before exponential/pole projection.
def raw(m,n):
 if m==n:return bd[m]-ad[m]+(4 if m==0 else 2)*mG[m]
 return 2*(d[m]-d[n])/(w[m]**2-w[n]**2)
rr=[(F(1,4)+w[p]**2)/(F(1,4)+w[n]**2) for n in ns[1:]]
R=[];ex=[(F(4*k+1,2),(-1)**(p+1)*2*F(4*k+1,2)**2*e2**(4*k+1)) for k in range(1,40)];t=[[1/(ak*ak+w[n]**2)-r/(ak*ak+w[p]**2) for n,r in zip(ns[1:],rr)] for ak,_ in ex];tail=(4*e2**161/(1-e2**4)).h
for m in ms:
 r0=raw(m,p);vals=[raw(m,n)-r*r0 for n,r in zip(ns[1:],rr)];
 for (ak,coef),v in zip(ex,t):
  x=coef/(ak*ak+w[m]**2);vals=[v0+x*v1 for v0,v1 in zip(vals,v)]
 R.append([v.widen(tail)/(I(2) if m==0 else sqrtint(2)) for v in vals])
R0,eR=grid(R); a=np.load(root/(par+'_normalization.npz')); B=obj(a['coefficient']);eB=F(str(a['coefficient_error']));Bscale=int(a['coefficient_scale']);Y=R0@B;YS=S*Bscale;eY=norm(R0,S)*eB+eR*(norm(B,Bscale)+eB)
# Reuse independently enclosed BN pilot overlaps and physical Gram.
bn=root.parent/'cr_int_t1d_bn'; old=np.load(bn/(par+'_head.npz')); assert np.array_equal(old['indices'],ms)
F0=obj(old['overlap']);FS=int(old['scale']);eF=F(str(old['overlap_error']));Gd=np.load(bk/par/'physical_gram.npz');G=obj(Gd['K']);GS=int(Gd['scale']);eG=F(str(Gd['error']))
# Trial physical projection onto F, only central head; residual norm bound using full physical G.
# y0 in orthonormal Neumann head; W approximates G^-1 F* y0. Exact norm squared central
# plus errors <= ||y0||_F^2 - 2 tr(W^* F0 Y0) + tr(W^* G0 W) + 2||W||_F eF ||y0||_F + eG||W||_F^2.
T=F0@Y;TS=FS*YS;V=rnd(np.linalg.solve(fl(G,GS),fl(T,TS)));nV=norm(V,S);nY=norm(Y,YS);trace=lambda a,b:sum(int(x)*int(y) for x,y in zip(a.flat,b.flat));YY=F(trace(Y,Y),YS**2);cross=F(trace(V,T),S*TS);vv=F(trace(V,G@V),GS*S*S);upper=YY-2*cross+vv+2*nV*eF*nY+eG*nV*nV
assert upper>=0
# bound square root to 1e-9 grid, exact.
def sqrtup(x):
 k=-((-x.numerator*10**100)//x.denominator);r=isqrt(k);return F(r+(r*r<k),10**50)
head=sqrtup(upper)+eY;rawbound=nY+eY
assert head<rawbound
print(par,'cutoff',cut,'rows',len(ms),'action enclosure error',float(eR),'physical overlap error',float(eF),'central head norm',float(nY),'central trial projection norm upper',float(sqrtup(upper)),'lift uncertainty allowance',float(eY),'projected head Frobenius upper',float(head),'raw bound',float(rawbound),flush=True)
np.savez_compressed(root/(par+'_head.npz'),cutoff=cut,indices=np.array(ms),action=np.array(R0,str),action_error=str(eR),overlap=np.array(F0,str),overlap_scale=str(FS),overlap_error=str(eF),scale=str(S),projection_trial=np.array(V,str),upper=str(head),rawupper=str(rawbound),central_norm=str(nY),central_projection_upper=str(sqrtup(upper)),lift_error=str(eY))
print(par,'BO energy-normalized projected finite-head pilot PASS (partial head only)',flush=True)
