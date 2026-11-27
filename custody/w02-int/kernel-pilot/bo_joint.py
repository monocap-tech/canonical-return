"""Enclosed exact KKT minimizer of signed energy with four BM moments, preserving BL constraints."""
from pathlib import Path
import sys,numpy as np
from fractions import Fraction as F
from math import isqrt
root=Path(__file__).resolve().parent;bl=root.parent/'cr_int_t1d_bl';bm=root.parent/'cr_int_t1d_bm';bn=root.parent/'cr_int_t1d_bn'
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2
p=int(sys.argv[1]);par=['even','odd'][p];S=10**80;obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1,upper=True):
 k=sum(int(v)**2 for v in a.flat);r=isqrt(k);return F(r+(upper and r*r<k),s)
def rnd(a):return np.vectorize(lambda v: int(F(float(v))*S),otypes=[object])(a)
def inv(A):
 n=len(A);I0=np.eye(n,dtype=object)*S; X=rnd(np.linalg.inv(fl(A)))
 for _ in range(5): X=X@(2*I0-A@X//S)//S
 return X
def fl(a,s=S):return np.array(a,float)/s
ci=lambda a,b:-((-a)//b)
def err(a,l,h,s):return norm(np.maximum(abs(a*D-l*s),abs(h*s-a*D)),s*D)
o=np.load(bl/(par+'_lift.npz'));n=np.load(root/(par+'_overlap.npz'));q=np.load(bn/(par+'_extension_energy.npz'));d=len(o['indices']);ns=np.r_[p,n['indices']];m=len(ns)-1-d
assert m==20 and np.array_equal(ns[1:],q['indices'])
pre=np.load(bn/(par+'_extension_precondition.npz'));assert F(str(pre['lower']))>F(999,1000)
s=int(n['scale']); HH=obj(n['H'])*S//s; Ho=HH[:,:d];Hn=HH[:,d:];eo=F(str(n['error']))+F(len(ns),S);en=eo
Z=inv(Ho);rho=norm(Z@Ho-np.eye(d,dtype=object)*S*S,S*S)+norm(Z,S)*eo;assert rho<F(1,100)
C=Z@Hn//S;nc=norm(C,S);ec=norm(Z,S)*(en+eo*nc+norm(Ho@C-Hn*S,S*S))/(1-rho)
print(par,'H inverse residual',float(rho),'C error',float(ec),flush=True)
lo=obj(q['lo']);hi=obj(q['hi']);Q=(lo+hi)*S//(2*D);A=Q[:d,:d];B=Q[:d,d:];T=Q[d:,d:];ea=err(A,lo[:d,:d],hi[:d,:d],S);eb=err(B,lo[:d,d:],hi[:d,d:],S);et=err(T,lo[d:,d:],hi[d:,d:],S);na=norm(A,S);nb=norm(B,S)
J=B.T*S-C.T@A;Du=T*S*S-(B.T@C+C.T@B)*S+C.T@A@C
# Round to common S grid and add exact truncation allowances.
J0=J//S;D0=Du//(S*S);ej=eb+ec*na+(nc+ec)*ea+norm(J-J0*S,S*S)
ed=et+2*(nb*ec+(nc+ec)*eb)+(2*nc*ec+ec*ec)*na+(nc+ec)**2*ea+norm(Du-D0*S*S,S**3)
W=4096;w=[pi*int(k)/4 for k in ns];assert all(x.h<W*D for x in w);pp={}
for a in range(2,55):
 if all(a%i for i in range(2,isqrt(a)+1)):
  k=a
  while k<=54:pp[k]=a;k*=a
wt=[(logi(I(a))/sqrtint(k),logi(I(k))) for k,a in pp.items()];assert len(wt)==24
h0=[];h1=[]
for idx,x in zip(ns,w):
 dj=x*psi(int(idx))[2]/2+sum((a*x*trig(x*s,1) for a,s in wt),I(0))
 h0.append(dj+sum(((-1)**p*F(4*k+1,2)**2*e2**(4*k+1)/(F(4*k+1,2)**2+x*x) for k in range(1,40)),I(0)))
 h1.append(dj*x*x-sum(((-1)**p*F(4*k+1,2)**4*e2**(4*k+1)/(F(4*k+1,2)**2+x*x) for k in range(1,40)),I(0)))
assert (e2**4*F(163,161)**2).h<D//2
h0=[v.widen((2*e2**161).h) for v in h0];h1=[v.widen((2*F(161,2)**2*e2**161).h) for v in h1]
rr=[(F(1,4)+w[0]**2)/(F(1,4)+x*x) for x in w[1:]];raw=[[I(1) for x in w],[x/W for x in h0],[x*x/W**2 for x in w],[x/W**3 for x in h1]]
mm=[[row[i+1]-rr[i]*row[0] for i in range(len(ns)-1)] for row in raw];ml=np.array([[v.l*S//D for v in row] for row in mm],object);mh=np.array([[ci(v.h*S,D) for v in row] for row in mm],object);M=(ml+mh)//2;em=norm(np.maximum(M-ml,mh-M),S)
Mo=M[:,:d];K=M[:,d:]*S-Mo@C;K0=K//S;ek=norm(K-K0*S,S*S)+norm(Mo,S)*ec+em*(1+nc+ec)
# saddle system D t + K^T lambda = -J; K t = -Mo. Scale moment multiplier by 2 (congruence).
scale=2;P=np.block([[D0,scale*K0.T],[scale*K0,np.zeros((4,4),dtype=object)]]);rhs=-np.vstack([J0,scale*Mo]);ep=ed+2*scale*ek;er=ej+scale*em
Y=inv(P);ny=norm(Y,S);res=norm(Y@P-np.eye(m+4,dtype=object)*S*S,S*S)+ny*ep
print(par,'C,J,D,K errors',*[float(v) for v in [ec,ej,ed,ek]],'KKT residual',float(res),'inverse norm',float(ny),flush=True)
assert res<F(1,2)
X=Y@rhs//S;nx=norm(X,S);delta=ny*(norm(P@X-rhs*S,S*S)+ep*nx+er)/(1-res)
# exact lift coefficients [(I-C t);t] and reference pole adjustment.
t=X[:m];nt=norm(t,S)+delta;rnorm=F(isqrt(sum(v.absmax()**2 for v in rr))+1,D);ceilroot=lambda n:isqrt(n)+(isqrt(n)**2<n)
a1=(ceilroot(d+m)+rnorm)*(1+(1+nc+ec)*nt)
# stronger coefficient bound directly with central corrected coefficient and propagated error.
B0=np.vstack([np.eye(d,dtype=object)*S*S-C@t,t*S]);eB=(1+nc+ec)*delta+ec*norm(t,S)
a1b=(ceilroot(d+m)+rnorm)*(norm(B0,S*S)+eB)
a1=min(a1,a1b)
old=F(str(np.load(bm/(par+'_moments.npz'))['a1bound']));assert a1<old
# exact true t gives unique constrained minimizer since D_U positive by inherited full finite positivity.
print(par,'KKT solution error',float(delta),'correction Frobenius bound',float(nt),'coefficient l1 bound',float(a1),'old',float(old),'ratio',float(a1/old),flush=True)
for N in [2**14,2**16,2**18,2**20,2**24]:
 tail=(74*W**4*a1)**2*F(4,3)**10/F(18*N**9);print('N',N,'tail squared',float(tail),flush=True)
np.savez_compressed(root/(par+'_joint.npz'),indices=ns[1:],scale=str(S),C=np.array(C,str),Cerror=str(ec),M=np.array(M,str),Merror=str(em),K=np.array(K0,str),Kerror=str(ek),D=np.array(D0,str),Derror=str(ed),J=np.array(J0,str),Jerror=str(ej),P=np.array(P,str),rhs=np.array(rhs,str),Y=np.array(Y,str),X=np.array(X,str),rho=str(res),delta=str(delta),a1bound=str(a1),coefficient=np.array(B0,str),coefficient_error=str(eB))
print(par,'BO exact constrained energy minimizer and improved tail PASS',flush=True)
