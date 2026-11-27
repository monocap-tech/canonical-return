"""High precision directed overlap evaluator; exact integer cross products."""
from pathlib import Path
import numpy as np, sys, hashlib
from fractions import Fraction as F
from math import isqrt
from bn_energy import I,D,pi,sqrtint
root=Path(__file__).resolve().parent; bk=root.parent/'cr_int_t1d_bk';bn=root.parent/'cr_int_t1d_bn'
p=int(sys.argv[1]);par=['even','odd'][p]; a=np.load(bk/par/'physical_dictionary.npz'); e=np.load(bk/par/'exceptional.npz');j=a['indices'];weights=a['weights'];E=e['E'].astype(object);ES=int(e['scale']);ws=int(a['weight_scale']);ns=np.load(bn/(par+'_indices.npy'));d=E.shape[1]
S=10**75; ci=lambda a,b:-((-a)//b)
# Directed grid at S; analytic sine table with exact signs.
const=[4/pi,2*sqrtint(2)/pi]; clo=[v.l*S//D for v in const];chi=[ci(v.h*S,D) for v in const]
def term(k):
 if k==0:return S,S
 t=k%8;sg=(1 if 1<=t<=3 else -1 if 5<=t<=7 else 0)*(1 if k>0 else -1);l,h=clo[t%2],chi[t%2];return (sg*l if sg>=0 else sg*h)//abs(k),ci(sg*h if sg>=0 else sg*l,abs(k))
def mul(l,h,a,b):
 v=(l*a,l*b,h*a,h*b);return min(v)//S,ci(max(v),S)
rt=isqrt(2*S*S);w=[isqrt(int(v)*S*S//ws) for v in weights]
def raw(jj,n,i):
 if not p and jj==0:l,h=(rt,rt+1) if n==0 else (0,0)
 else:
  l,h=term(jj-2*n);u,v=term(jj+2*n);l,h=(l+u,h+v) if not p else (l-v,h-u)
  if ((n+p)//2)%2:l,h=-h,-l
 return mul(l,h,w[i],w[i]+1)
rr=[]
for n in ns:
 v=(F(1,4)+(p*pi/4)**2)/(F(1,4)+(int(n)*pi/4)**2);rr.append((v.l*S//D,ci(v.h*S,D)))
mid=np.empty((len(j),len(ns)),object);maxrad=0
for i,jj in enumerate(j):
 l0,h0=raw(int(jj),p,i)
 for k,n in enumerate(ns):
  l,h=raw(int(jj),int(n),i);a,b=mul(l0,h0,*rr[k]);l,h=l-b,h-a;m=(l+h)//2;mid[i,k]=m;maxrad=max(maxrad,m-l,h-m)
H=E.T@mid; h0=H//ES
# absolute row coefficient sum bounds entry enclosure; Frobenius norm bound via max entry.
maxsum=max(sum(abs(int(v)) for v in E[:,k]) for k in range(d));er=F(maxsum*maxrad,ES*S)+F(1,S);eta=max(d,len(ns))*er
# Check compatibility with the previous independently enclosed centers.
o=np.load(root.parent/'cr_int_t1d_bl'/(par+'_lift.npz')); n=np.load(bn/(par+'_extension_lift.npz'));old=np.concatenate([o['H'],n['Hnew']],axis=1);ss=int(o['scale']);eo=max(F(str(o['error'])),F(str(n['error'])));assert all(abs(F(int(x),S)-F(int(y),ss))<=eo+er for x,y in zip(h0.flat,old.flat))
np.savez_compressed(root/(par+'_overlap.npz'),H=np.array(h0,str),scale=str(S),error=str(eta),indices=ns,entry_error=str(er));print(par,'entries',h0.size,'raw interval radius',maxrad,'cross error',float(eta),'old compatibility PASS',flush=True)
