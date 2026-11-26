from pathlib import Path
from fractions import Fraction as F
from math import isqrt,factorial
import sys,numpy as np,time
root=Path(__file__).resolve().parent;parent=root.parent;bp=parent/'cr_int_t1d_bp';S=10**12
src=(bp/'bp_project.py').read_text();exec(src[src.index('obj=lambda'):src.index('h=np.load(root/')]);a=np.load(root/'bs_trial.npz');assert int(a['scale'])==S
# all finite trigonometric terms expressed as coefficients of unweighted sine(x*pi*j/8)/2 on I
js=a['js'];ms=a['ms'];yn=a['yn'];gc=a['gc'];c=a['c'];r=isqrt(2*S*S);sign=np.where(((ms+1)//2)%2,-1,1);cy=(yn.astype(object)*r*sign[:,None]//S).astype(np.int64)
jall=np.r_[js,2*ms];C=np.vstack([c,cy]);assert len(set(jall))==len(jall) and np.all(jall>0)
# additional conversion error sqrt(2) plus floor; use contraction separately on odd/even sets
err=np.array([F(str(e))+norm(yn[:,k:k+1],S)/S+F(len(ms),S) for k,e in enumerate(a['errors'])],object)
# Exact rational pi enclosure from Machin; sine residues exact via sqrt2.
atan=lambda x,k:sum(((-1)**i*x**(2*i+1)/F(2*i+1) for i in range(k)),F(0))
pl=16*atan(F(1,5),40)-4*atan(F(1,239),39);ph=16*atan(F(1,5),39)-4*atan(F(1,239),40);assert 3<pl<ph<F(22,7)
lo=[int(2*S/ph),int(F(r)/ph)];hi=[-((-F(2*S,1)/pl)//1),-((-F(r+1,1)/pl)//1)];assert all(0<b-a0<=2 for a0,b in zip(lo,hi));lo=np.array(lo,np.int64);hi=np.array(hi,np.int64)
def term(k):
 m=k%8;sg=np.where((m>=1)&(m<=3),1,np.where((m>=5)&(m<=7),-1,0))*np.where(k<0,-1,1);ix=m%2;den=np.maximum(abs(k),1);l=np.where(sg>=0,sg*lo[ix],sg*hi[ix])//den;h=-((-np.where(sg>=0,sg*hi[ix],sg*lo[ix]))//den);return l,h
# g_j = <u_j,sinh(x/2)> from BJ primitive. Use high precision rational interval from Taylor.
e=sum((F(1,factorial(i)) for i in range(80)),F(0));eh=e+F(2,factorial(80));shl=(e-1/e)/2;shh=(eh-1/eh)/2;chl=(e+1/e)/2;chh=(eh+1/eh)/2
rt=(F(r,2*S),F(r+1,2*S));z=(F(0),F(0));one=(F(1),F(1));neg=lambda a:(-a[1],-a[0]);cs=[one,rt,z,neg(rt),neg(one),neg(rt),z,rt];sn=[z,rt,one,rt,z,neg(rt),neg(one),neg(rt)]
def mul(a,b):v=[x*y for x in a for y in b];return min(v),max(v)
def gp(j):
 w=(pl*int(j)/8,ph*int(j)/8);c0=cs[int(j)%8];s0=sn[int(j)%8];u=mul((chl/2,chh/2),s0);v=mul(mul(w,(shl,shh)),c0);num=(u[0]-v[1],u[1]-v[0]);den=(F(1,4)+w[0]**2,F(1,4)+w[1]**2);l,h=mul(num,(1/den[1],1/den[0]));return int(l*S//1),int(-((-h*S)//1))
# precompute g overlaps using a 1e18 directed integer interval implementation to avoid large fractions per j
# equivalent to primitive; integer floor/ceil every op.
D=10**18
def fi(a):return int(a[0]*D//1),int(-((-a[1]*D)//1))
def im(a,b):v=[i*j for i in a for j in b];return min(v)//D,-((-max(v))//D)
def ig(j):
 w=fi((pl*int(j)/8,ph*int(j)/8));u=im(fi((chl/2,chh/2)),fi(sn[int(j)%8]));v=im(im(w,fi((shl,shh))),fi(cs[int(j)%8]));num=(u[0]-v[1],u[1]-v[0]);den=im(w,w);den=(den[0]+D//4,den[1]+D//4);assert den[0]>0
 vals=[F(i*S,k) for i in num for k in den];return int(min(vals)//1),int(-((-max(vals))//1))
# cutoff is a common circular frequency index, not physical Neumann index.
cut=int(sys.argv[1]) if len(sys.argv)>1 else 32768;ids=np.union1d(jall,np.arange(1,cut,2,dtype=np.int64));mode=sys.argv[2] if len(sys.argv)>2 else 'all';start=time.time();print('circle output rows',len(ids),'inputs',len(jall),flush=True)
# directed exact Gram action; each entry grid radius <=4; use exact contraction of center, and per-column l1 for a row error
# use coefficient grid fixed S, g grid S; output grid S^2. g term explicit.
l1=np.array([sum(abs(int(x)) for x in col) for col in C.T],object);Gout=[];erout=[];gr=[]
part=0;startrow=0;endrow=0;
for b in range(startrow,endrow,256):
 j=ids[b:min(b+256,endrow)];d=j[:,None]-jall[None,:];s=j[:,None]+jall[None,:];l,h=term(d);u,v=term(s);l=l-v+(d==0)*(S//2);h=h-u+(d==0)*(S//2);mid=(l+h)//2;rad=np.maximum(mid-l,h-mid);assert rad.max()<=4
 # int64 input checked by dot helper; no approximate cancellation assumptions
 Z=dot(mid,C);eZ=dot(rad,abs(C));gints=np.array([ig(int(x)) for x in j],np.int64);gm=gints.sum(1)//2;ger=np.maximum(gm-gints[:,0],gints[:,1]-gm);Z+=gm[:,None].astype(object)@gc.astype(object);eZ+=ger[:,None].astype(object)@abs(gc).astype(object)
 Gout.append(Z);erout.append(eZ);gr.extend(gints.tolist())
 if b%2048==0:print('through',int(j[-1]),'elapsed',round(time.time()-start,1),flush=True)
parts=[np.load(root/f'circle_part{k}.npz') for k in range((len(ids)+8191)//8192)];assert np.array_equal(np.concatenate([p['ids'] for p in parts]),ids);Z=np.concatenate([obj(p['Z']) for p in parts]);eZ=np.concatenate([obj(p['Zerr']) for p in parts]);gr=np.concatenate([p['gr'] for p in parts]);gp0=gr.sum(1)//2;gper=np.maximum(gp0-gr[:,0],gr[:,1]-gp0);nn=np.searchsorted(ids,jall);assert np.array_equal(ids[nn],jall)
# v0 norm gram: C*G(C+g) + gc*(g*C+||g||² gc)
G1=C.astype(object).T@Z[nn];G2=gc.astype(object).T@(gp0[nn][None,:].astype(object)@C.astype(object));bk=parent/'cr_int_t1d_bk'/'odd';pd=np.load(bk/'pole.npz');gl=F(str(pd['gnormlo']));gh=F(str(pd['gnormhi']));gmid=int((gl+gh)*S/2);egn=max(abs(F(gmid,S)-gl),abs(F(gmid,S)-gh));G=G1+G2+gmid*(gc.astype(object).T@gc.astype(object));G=(G+G.T)//2;# scale S^3
# entrywise error from Gram action and g vector, conservative symmetrization; store per-entry matrix
E1=abs(C).astype(object).T@eZ[nn];E2=abs(gc).astype(object).T@(gper[nn][None,:].astype(object)@abs(C).astype(object));E3=abs(gc).astype(object).T@abs(gc).astype(object)*int(-((-egn*S)//1));Er=np.maximum(E1+E2+E3,(E1+E2+E3).T)+1
jj=np.flatnonzero((ids%2==1)&(ids<cut));assert np.array_equal(ids[jj],np.arange(1,cut,2));
np.savez_compressed(root/'bs_circle.npz',ids=ids[jj],Z=np.array(Z[jj],str),Zerr=np.array(eZ[jj],str),Zscale=str(S*S),G=np.array(G,str),Gerr=np.array(Er,str),Gscale=str(S**3),reprerr=np.array([str(e) for e in err]),cutoff=cut)
print('circle exact Gram and partial action PASS',len(jj),'elapsed',round(time.time()-start,1),flush=True)
