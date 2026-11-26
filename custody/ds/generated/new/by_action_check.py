"""Independent recorded-proposal recheck, inherited custody and signed row checks."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,numpy as np
from bn_energy import I,D,pi,psi,logi,sqrtint,trig,e2
root=Path(__file__).resolve().parent;bp=root.parent/'cr_int_t1d_bp';bo=root.parent/'cr_int_t1d_bo';bk=root.parent/'cr_int_t1d_bk';bn=root.parent/'cr_int_t1d_bn';p=1;par='odd'
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 k=sum(int(v)**2 for v in a.flat);r=isqrt(k);return F(r+(r*r<k),s)
def sqrtup(x,s=10**40):
 assert x>=0;k=-((-x.numerator*s*s)//x.denominator);r=isqrt(k);return F(r+(r*r<k),s)
def rows(a,s=1):return F(max(sum(abs(int(v)) for v in row) for row in a),s)
def lower(a,s=1):return F(min(int(a[i,i])-sum(abs(int(a[i,j])) for j in range(len(a)) if j!=i) for i in range(len(a))),s)
a=np.load(root/'odd_corrector.npz');j=np.load(bo/'odd_joint.npz');ns=[1]+list(map(int,a['indices']));assert np.array_equal(a['indices'][:len(j['indices'])],j['indices']);assert len(set(ns))==len(ns);assert len(a['new'])==16
S=int(a['coefficient_scale']);Z=obj(a['Z']);zs=int(a['Zscale']);P=obj(a['P']);Pn=obj(a['Pn']);ep=F(str(a['Perror']));en=F(str(a['Pnerror']));C=obj(a['C']);U=obj(a['coefficient']);nz=norm(Z,zs);rho=norm(P@Z-np.eye(len(P),dtype=object)*S*zs,S*zs)+ep*nz;assert rho==F(str(a['rho']))<F(1,100);Cbar=Z@Pn;rc=norm(Cbar-C*zs,zs*S);assert rc==F(str(a['roundC']));err=nz*(norm(P@Cbar-Pn*S*zs,S*S*zs)+ep*norm(Cbar,zs*S)+en)/(1-rho)+rc;assert err==F(str(a['coefficient_error']))<F(1,10**30);assert np.array_equal(U,np.vstack([-C,np.eye(16,dtype=object)*S]));print('exact range-restricted right inverse and center rounding PASS',flush=True)
# New full action coverage and input-allowance compatibility with the executed chunk runs.
h=np.load(root/'odd_U_full.npz');hw=np.load(bp/'odd_head_full.npz');ms=h['indices'];assert np.array_equal(ms,hw['indices']) and np.array_equal(ms,np.arange(1,16385,2));Y=obj(h['Y']);YS=int(h['scale']);ey=F(str(h['error']));assert ey==sum(map(F,map(str,h['tile_errors'])),F(0))<F(1,10**5)*4;chunks=[np.load(root/f'odd_U_{lo}_{hi}.npz') for lo,hi in [(1,2047),(2049,6143),(6145,11263),(11265,16383)]];assert np.array_equal(h['Y'],np.concatenate([c['Y'] for c in chunks]));assert np.array_equal(ms,np.concatenate([c['indices'] for c in chunks]));assert ey==sum((F(str(c['error'])) for c in chunks),F(0));assert err+F(U.size,10**65)<F(1,10**45);print('full head coverage and chunk input allowance PASS',flush=True)
# Independently evaluate old-style signed divided difference and individual exponential terms
# on selected rows, including reference and new trial modes.
pp={}
for q in range(2,55):
 if all(q%i for i in range(2,isqrt(q)+1)):
  k=q
  while k<=54:pp[k]=q;k*=q
assert len(pp)==24
wt=[(logi(I(q))/sqrtint(k),logi(I(k))) for k,q in pp.items()];test=[1,257,773,2039,8193,16383];alln=sorted(set(ns+test));w={n:psi(n)[0] for n in alln};J={n:w[n]*psi(n)[2]/2 for n in alln};bd={n:psi(n)[2]/(2*w[n])+psi(n)[3]/4 for n in alln};mG={n:psi(n)[1]-logi(pi) for n in alln};sn=[{n:trig(w[n]*s,1) for n in alln} for _,s in wt];cs=[{n:trig(w[n]*s,0) for n in alln} for _,s in wt];ad={n:sum((a*((4-s)*c[n]-v[n]/w[n]) for (a,s),v,c in zip(wt,sn,cs)),I(0)) for n in alln};d={n:J[n]+w[n]*sum((a*v[n] for (a,_),v in zip(wt,sn)),I(0)) for n in alln}
def raw(m,n):
 if m==n:return bd[m]-ad[m]+2*mG[m]
 return 2*(d[m]-d[n])/(w[m]**2-w[n]**2)
rr=[(F(1,4)+w[p]**2)/(F(1,4)+w[n]**2) for n in ns[1:]];ex=[(F(4*k+1,2),2*F(4*k+1,2)**2*e2**(4*k+1)) for k in range(1,40)];tr=[[1/(ak*ak+w[n]**2)-r/(ak*ak+w[p]**2) for n,r in zip(ns[1:],rr)] for ak,_ in ex];tail=(4*e2**161/(1-e2**4)).h;S2=10**75
for m in test:
 vals=[raw(m,n)-r*raw(m,p) for n,r in zip(ns[1:],rr)]
 for (ak,c),row in zip(ex,tr):
  x=c/(ak*ak+w[m]**2);vals=[v+x*t for v,t in zip(vals,row)]
 vals=[v.widen(tail)/sqrtint(2) for v in vals];lo=np.array([v.l*S2//D for v in vals],object);hi=np.array([-((-v.h*S2)//D) for v in vals],object);R=(lo+hi)//2;er=norm(np.maximum(R-lo,hi-R),S2);y=R@U;eta=norm(R,S2)*err+er*(norm(U,S)+err);k=int(np.where(ms==m)[0][0]);assert norm(y*YS-Y[k]*S2*S,S2*S*YS)<=eta+ey,m
print(len(test)*U.shape[1],'independent direct signed action entry overlaps PASS',flush=True)
