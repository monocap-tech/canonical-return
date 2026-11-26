"""Independent exact check of CL's stored right-inverse certificate.
No historical numerical threshold is assumed, except the input allowance
actually used by the inherited directed action and the source construction."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np, sys, json
a=np.load(sys.argv[1]);j=np.load(sys.argv[2]);old=np.load(sys.argv[3])
obj=lambda x:np.vectorize(int,otypes=[object])(x)
def norm(x,s=1):
    q=sum(int(v)**2 for v in x.flat);r=isqrt(q);return F(r+(r*r<q),s)
S=int(a['coefficient_scale']);zs=int(a['Zscale'])
assert S==10**60 and zs==S*S
idx=a['indices'];new=a['new'];assert len(new)==len(set(new))==16
assert np.array_equal(idx,np.r_[j['indices'],new])
assert len(set(idx))==len(idx) and all(idx%2==1) and 1 not in idx
assert not np.intersect1d(new,old['indices']).size
P,Pn,Z,C,U=map(obj,[a['P'],a['Pn'],a['Z'],a['C'],a['coefficient']])
ep,en=F(str(a['Perror'])),F(str(a['Pnerror']))
assert ep>=0 and en>=0 and Pn.shape==(P.shape[0],16)
assert Z.shape==(P.shape[1],P.shape[0]) and C.shape==(P.shape[1],16)
nz=norm(Z,zs)
rho=norm(P@Z-np.eye(len(P),dtype=object)*S*zs,S*zs)+ep*nz
assert rho==F(str(a['rho']))<F(1,100)
Cbar=Z@Pn;rc=norm(Cbar-C*zs,zs*S)
assert rc==F(str(a['roundC']))
err=nz*(norm(P@Cbar-Pn*S*zs,S*S*zs)+ep*norm(Cbar,zs*S)+en)/(1-rho)+rc
assert err==F(str(a['coefficient_error']))==F(str(a['error']))<F(1,10**30)
assert np.array_equal(U,np.vstack([-C,np.eye(16,dtype=object)*S]))
# bq_action computes its own outward-rounded input allowance, not a
# fixed historical coefficient-error threshold. Reproduce that allowance.
BS=10**65;grid=10**45
allow=err+F(U.size,BS);eB=F(-((-allow.numerator*grid)//allow.denominator),grid)
assert eB>=allow
rec=dict(columns=16,rank_receiver=P.shape[0],rho=str(rho),coefficient_error=str(err),action_input_allowance=str(eB),new_indices=list(map(int,new)))
print('CM EXACT RIGHT INVERSE AUDIT',json.dumps(rec),flush=True)
if len(sys.argv)>4:Path(sys.argv[4]).write_text(json.dumps(rec,indent=2)+'\n')
