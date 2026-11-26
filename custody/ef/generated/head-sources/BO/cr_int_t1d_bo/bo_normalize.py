"""Congruence certificate for the true constrained-lift signed energy.
Exact integer center products, inherited interval matrix and Frobenius error."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,sys
root=Path(__file__).resolve().parent;bn=root.parent/'cr_int_t1d_bn';p=int(sys.argv[1]);par=['even','odd'][p]
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 t=sum(int(v)**2 for v in a.flat);r=isqrt(t);return F(r+(r*r<t),s)
a=np.load(root/(par+'_joint.npz'));q=np.load(bn/(par+'_extension_energy.npz'));s=int(a['scale'])**2;B=obj(a['coefficient']);eb=F(str(a['coefficient_error']));D=int(q['scale']);lo=obj(q['lo']);hi=obj(q['hi']);M=(lo+hi)//2;rad=np.maximum(M-lo,hi-M);n=B.shape[1];S=10**110
# central energy, round symmetric matrix to S grid.
T=B.T@M@B; den=s*s*D; A=T*S//den; eta=norm(T*S-A*den,S*den)+norm(M,D)*(2*norm(B,s)*eb+eb*eb)+norm(rad,D)*(norm(B,s)+eb)**2
# exact fixed point LDL proposal, then verified congruence (not relying on LDL accuracy)
P=A.copy();L=np.eye(n,dtype=object)*S;pv=[]
for i in range(n):
 v=int(P[i,i]);pv.append(v);assert v>0
 for j in range(i+1,n):L[j,i]=int(P[j,i])*S//v
 for j in range(i+1,n):
  for k in range(j,n):P[k,j]=P[j,k]=int(P[k,j])-int(P[k,i])*int(L[j,i])//S
J=np.eye(n,dtype=object)*S
for i in range(n):
 for j in range(i):J[i,j]=-sum(int(L[i,k])*int(J[k,j]) for k in range(j,i))//S
Z=J*np.array([isqrt(S**3//v) for v in pv],object)[:,None]//S
assert all(Z[i,i]>0 for i in range(n)) and all(Z[i,j]==0 for i in range(n) for j in range(i+1,n))
Q=Z@A@Z.T; nZ=norm(Z,S);e=eta*nZ*nZ
row=[sum(abs(int(Q[i,j])-(S**3 if i==j else 0)) for j in range(n)) for i in range(n)];r=F(max(row),S**3)+e
print(par,'coefficient error',float(eb),'energy error',float(eta),'min midpoint pivot',float(F(min(pv),S)),'normalizer norm',float(nZ),'normalized energy error',float(r),flush=True)
assert r<F(1,100)
# exact centered normalized coefficient product, with true norm error.
C=B@Z.T; CS=s*S;ec=eb*nZ;print(par,'normalized coefficient norm',float(norm(C,CS)),'normalized coefficient error',float(ec),'PASS',flush=True)
np.savez_compressed(root/(par+'_normalization.npz'),Z=np.array(Z,str),scale=str(S),energy_error=str(eta),relative_error=str(r),coefficient=np.array(C,str),coefficient_scale=str(CS),coefficient_error=str(ec))
