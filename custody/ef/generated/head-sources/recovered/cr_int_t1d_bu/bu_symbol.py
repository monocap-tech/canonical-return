from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,sys
root=Path(__file__).resolve().parent;bi=root.parent/'cr_int_t1d_bi'/'even';sys.path.insert(0,str(bi));import os;os.chdir(bi);import ax_arithmetic as ar
j=np.load(root/'surplus_circle.npz')['ids'];assert np.array_equal(j,np.arange(0,32768,2));G=10**9;al=np.zeros(len(j),np.int64);au=al.copy();sumupper=F(0)
for n,base in ar.p.pp.items():
 l,h=ar.p.logn(base);r=isqrt(n*10**40);wl=int(l*10**20*G//(r+1));wh=int(-((-h*10**20*G)//r));assert wh-wl<=2 and 0<wl<=wh<G;sumupper+=F(2*wh,G)
 c,s=ar.phases(n,j);a,b=ar.mul(c-4,c+4,np.array(wl),np.array(wh));al+=2*a;au+=2*b
assert len(ar.p.pp)==24 and sumupper<25
bmin=F(6,5)*(F(399,50)-sumupper);assert bmin<0
# Lower nonnegative symbol b_j-b_min=(6/5)(sumupper-a_j)
w=np.maximum(0,(6*(int(sumupper*G)-au))//5);assert w.max()<60*G and np.all(al<=au)
np.savez_compressed(root/'bu_symbol.npz',ids=j,w=w,scale=G,bmin=str(bmin),sumupper=str(sumupper),al=al,au=au)
print('signed symbol PASS',len(j),'sum upper',sumupper,'bmin',bmin,float(bmin),flush=True)
