"""Full head physical projection with directed overlaps and radix-split exact dots.
All float GEMMs only operate on integer limbs with checked absolute dot bounds <2^53,
under the inherited ordinary binary64 arithmetic model. Float solve is only a trial.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,factorial
import sys,numpy as np,time
from bn_energy import I,D as ID,pi,sqrtint
root=Path(__file__).resolve().parent;bk=root.parent/'cr_int_t1d_bk';bn=root.parent/'cr_int_t1d_bn';bo=root.parent/'cr_int_t1d_bo';p=int(sys.argv[1]);par=['even','odd'][p];S=10**12
obj=lambda a:np.vectorize(int,otypes=[object])(a)
def norm(a,s=1):
 k=sum(int(v)**2 for v in a.flat);r=isqrt(k);return F(r+(r*r<k),s)
def up(x,s=10**40):return F(-((-x.numerator*s)//x.denominator),s)
def sqrtup(x,s=10**30):
 assert x>=0;k=-((-x.numerator*s*s)//x.denominator);r=isqrt(k);return F(r+(r*r<k),s)
def fl(a,s):return np.vectorize(lambda v:float(F(int(v),s)),otypes=[float])(a)
def dot(A,B):
 # signed limbs, exact ordinary binary64 products and sums, with no cancellation assumptions.
 assert A.shape[1]==B.shape[0];r=2**16;assert A.shape[1]*(r-1)**2<2**53
 def limbs(a):
  a=np.asarray(a,dtype=np.int64);sg=np.sign(a);v=abs(a);out=[]
  while np.any(v):out.append((sg*(v%r)).astype(float));v//=r
  return out or [np.zeros(a.shape)]
 a=limbs(A);b=limbs(B);z=np.zeros((A.shape[0],B.shape[1]),object)
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   t=x@y;assert np.all(t==np.trunc(t));z+=t.astype(np.int64).astype(object)*r**(i+j)
 return z
bp=root.parent/'cr_int_t1d_bp';h=np.load(bp/(par+'_head_full.npz'));u=np.load(root/(par+'_U_full.npz'));assert np.array_equal(h['indices'],u['indices']) and int(h['scale'])==int(u['scale']);Y=np.concatenate([h['Y'],u['Y']],axis=1);YS=int(h['scale']);eY=F(str(h['error']));eU=F(str(u['error']));ms=h['indices'];nY=norm(Y,YS);d0=h['Y'].shape[1];assert ms[-1]>=16383
physical=np.load(bk/par/'physical_dictionary.npz');excp=np.load(bk/par/'exceptional.npz');pd=np.load(bk/par/'pole.npz');gd=np.load(bk/par/'physical_gram.npz');j=physical['indices'];weights=physical['weights'];E=excp['E'].astype(np.int64);ES=int(excp['scale']);pg=obj(gd['pole_compressed']);GS=int(gd['scale']);G=obj(gd['K']);eG=F(str(gd['error']));ep=F(str(gd['pole_error']));gl=F(str(pd['gnormlo']));gh=F(str(pd['gnormhi']))
# Inherited directed int64 overlap grid (D=1e9); exact products split into limbs.
src=(bn/'bn_lifts.py').read_text();s0=src[src.index('# directed overlaps'):src.index('# raw overlaps are weighted')];env={'np':np,'F':F,'isqrt':isqrt,'j':j,'weights':weights,'p':p};exec(s0,env);Dov=env['D']
e=I(sum((F(1,factorial(k)) for k in range(100)),F(0))).widen(I(F(2,factorial(100))).h);npg=norm(pg,GS)
Ftiles=[];errs=[];start=time.time();old=np.load(bo/(par+'_head.npz'));oldF=obj(old['overlap']);oldFS=int(old['overlap_scale']);oldErr=F(str(old['overlap_error']));nold=len(old['indices']);Ftiles.append(((oldF*S)//oldFS).astype(np.int64));errs.append(oldErr+F(oldF.size,S))
for begin in range(nold,len(ms),256):
 tile=ms[begin:begin+256];l,h=env['raw'](tile);mid=(l+h)//2;rad=np.maximum(mid-l,h-mid);H=dot(E.T,mid);Hr=dot(abs(E).T,rad);HS=Dov*ES
 # exact integer raw overlap projection into Neumann orthonormal coordinates
 mass=sqrtint(2);r=I(1)/mass;rm=(r.l+r.h)//2;er=F(max(rm-r.l,r.h-rm),ID)
 A0=(H*rm*S)//(HS*ID);eA=norm(Hr,HS)*(F(rm,ID)+er)+er*norm(H,HS)+F(A0.size,S)
 # subtract pole projection (same formula as BO); compute central g/gnorm row in intervals
 gr=[(((-1)**p*e-1/e)/2)/(F(1,4)+(pi*int(m)/4)**2)/mass/I(gl,gh) for m in tile]
 gm=np.array([(v.l*S//ID+(-((-v.h*S)//ID)))//2 for v in gr],object);ge=norm(np.array([max(abs(int(a)*ID-v.l*S),abs(int(a)*ID-v.h*S)) for a,v in zip(gm,gr)],object),S*ID);ng=norm(gm,S)
 pole=(pg[:,None]*gm[None,:])//GS;ePole=npg*ge+ep*(ng+ge)+F(pole.size,S)
 F0=A0-pole;assert max(abs(v) for v in F0.flat)<2**63;Ftiles.append(F0.astype(np.int64));errs.append(up(eA+ePole))
 if (begin-nold)%1024==0:print(par,'projection overlap through',tile[-1],'error',float(errs[-1]),'elapsed',round(time.time()-start,1),flush=True)
F0=np.concatenate(Ftiles,axis=1);assert F0.shape[1]==len(ms);eF=sum(errs,F(0));print(par,'full overlap error',float(eF),flush=True)
# Full physical projection: arbitrary F and Mellin-pole trials. F*g=0 exactly.
T=dot(F0,Y);TS=S*YS;V=np.rint(np.linalg.solve(fl(G,GS),fl(T,TS))*S).astype(np.int64)
# Actual orthonormal Neumann pairings with g. The n=0 mode has mass 4.
gr=[(((-1)**p*e-1/e)/2)/(F(1,4)+(pi*int(m)/4)**2)/(I(2) if m==0 else sqrtint(2)) for m in ms]
g0=np.array([(v.l*S//ID+(-((-v.h*S)//ID)))//2 for v in gr],object);eg=norm(np.array([max(abs(int(a)*ID-v.l*S),abs(int(a)*ID-v.h*S)) for a,v in zip(g0,gr)],object),S*ID)
gn=(gl+gh)/2;gn0=int(gn*GS);egn=max(abs(F(gn0,GS)-gl),abs(F(gn0,GS)-gh));assert gl>0
Tg=g0[None,:]@Y.astype(object);v=(Tg*GS)//(gn0*S);assert max(abs(int(x)) for x in v.flat)<2**63
YY=dot(Y.T,Y);Vo=V.astype(object);C=Vo.T@T;VV=Vo.T@G@Vo;Cg=v.T@Tg;VVg=gn0*(v.T@v);assert YS==S
P=YY*GS-C*GS//S-C.T*GS//S+VV-Cg*GS//S-Cg.T*GS//S+VVg;PS=S*S*GS;assert np.array_equal(P,P.T)
print(par,'combined physical projection PASS',P.shape,'eF',float(eF),'eg',float(eg),'egn',float(egn),flush=True)
np.savez_compressed(root/(par+'_combined_projection.npz'),P=np.array(P,str),scale=str(PS),F_error=str(eF),G_error=str(eG),g_error=str(eg),gn_error=str(egn),trial=V,pole_trial=np.array(v,str),trial_scale=str(S),Y=Y,Yscale=str(YS),head_error=str(eY),U_error=str(eU),d0=d0)
