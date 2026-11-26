from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np, hashlib, sys
root=Path(__file__).resolve().parent;parent=root.parent;bp=parent/'cr_int_t1d_bp';S=10**12
src=(bp/'bp_project.py').read_text();exec(src[src.index('obj=lambda'):src.index('h=np.load(root/')]);src=(parent/'cr_int_t1d_bt'/'bt_physical.py').read_text();exec(src[src.index('def rows('):src.index('bs=parent/')]);
# matrix interval helpers
sym=lambda a:(a+a.T)/2
I=lambda n: np.eye(n,dtype=object)
def cf(a,s):return np.vectorize(lambda x:F(int(x),s),otypes=[object])(a)
def scalarfloor(x,s=S):return int(x*s//1)
def scalarceil(x,s=S):return int(-((-x*s)//1))
def gram(a,e,w=None):
 if w is None:w=np.ones(a.shape[0],dtype=object)
 return a.T@(w[:,None]*a),abs(a).T@(abs(w)[:,None]*e)+e.T@(abs(w)[:,None]*abs(a))+e.T@(abs(w)[:,None]*e)
def load(par,tag):
 if par=='even':
  d=parent/'cr_int_t1d_bu';a=np.load(d/'bu_trial.npz');c=np.load(d/f'{tag}_circle.npz');o=np.load(d/f'{tag}_outside.npz');b=np.load(d/'bu_symbol.npz');en=np.load(d/'even_BQ_energy.npz');q=obj(a['Q']);X=obj(a['X']);x=obj(a['x']);t=obj(c['t']);alpha=obj(c['alpha']);eY=F(str(a['head_error']));eU=F(str(a['Uerror']));tau=F(str(a['tail']));tauU=F(str(a['Utail']));mu=F(str(a['mu']));eps=F(str(a['eps']));sl=slice(None)
 else:
  d=parent/'cr_int_t1d_bs';bq=parent/'cr_int_t1d_bq';a=np.load(d/'bs_trial.npz');c=np.load(d/'bs_circle.npz');o=np.load(parent/'cr_int_t1d_bt'/'bt_outside.npz');b=np.load(d/'bs_symbol.npz');en=np.load(bq/'odd_BQ_energy.npz');sl=slice(4,12);q=obj(a['Q'])[:,sl];X=obj(a['X']);x=X@q;t=obj(a['t'])[:,sl];alpha=obj(a['a']);A=np.load(bq/'odd_combined_projection.npz');B=np.load(bp/'odd_projected.npz');eY=F(str(A['head_error']));eU=F(str(A['U_error']));tau=F(str(B['tail']));T=np.load(bq/'odd_U_tail.npz');tauU=F(str(T['tail'][list(T['cutoff']).index(16384)]));mu=F(str(B['mu']));eps=F(str(B['eps']))
 K=obj(en['K']);D=obj(en['D']);KS=int(en['scale']);eK=F(str(en['Kerror']));eD=F(str(en['Derror']));nq=sqrtup(rows(cf(q.T@q,S*S)));nx=norm(x,S*S);E=(1-eps)*cf(q.T@q,S*S)+cf((x.T@K@q+q.T@K.T@x)*S+x.T@D@x,KS*S**4);eE=2*nx*eK*nq+nx*nx*eD;derr=(eY+sqrtup(tau))*nq+(eU+sqrtup(tauU))*nx
 diag=cf(t.T@((10**8-alpha)[:,None]*t),10**8*S*S)
 Z=obj(c['Z'])[:,sl];Ze=obj(c['Zerr'])[:,sl];ZS=int(c['Zscale']);G=cf(obj(c['G'])[sl,sl],int(c['Gscale']));Ge=cf(obj(c['Gerr'])[sl,sl],int(c['Gscale']));er=[F(str(e)) for e in c['reprerr'][sl]];e=sqrtup(sum(v*v for v in er));assert np.array_equal(c['ids'],b['ids']);WS=int(b['scale']);bl=[F(6,5)*(F(399,50)-F(int(v),WS)) for v in b['au']];bh=[F(6,5)*(F(399,50)-F(int(v),WS)) for v in b['al']];bmin=F(str(b['bmin']));bmax=F(6,5)*(F(399,50)+F(str(b['sumupper'])));bn=F(2394,125);assert all(l<=h for l,h in zip(bl,bh));bmin=min(bmin,min(bl));bmax=max(bmax,max(bh));
 # actual B gram center, with global signed-symbol tail interval rather than assume sign. Tail is G-Z*Z.
 zz,zze=gram(Z,Ze);Z2=cf(zz,ZS**2);Z2E=cf(zze,ZS**2);Vtail=G-Z2;et=rows(Ge+Z2E);vt=rows(Vtail)+et;assert vt>=0
 w=np.array([scalarfloor((l+h)/2) for l,h in zip(bl,bh)],object);we=max(max(abs(F(int(v),S)-l),abs(F(int(v),S)-h)) for v,l,h in zip(w,bl,bh));g,ge=gram(Z,Ze,w);mid=(bmin+bmax)/2;rad=(bmax-bmin)/2;Bg=cf(g,S*ZS**2)+mid*Vtail;Be=rows(cf(ge,S*ZS**2))+we*(rows(Z2)+rows(Z2E))+abs(mid)*et+rad*vt
 # B² upper using circular J² and physical exterior correction. Keep central tail matrix (not scalar).
 w2=np.array([scalarceil(max(l*l,h*h)) for l,h in zip(bl,bh)],object);g,ge=gram(Z,Ze,w2);J2=cf(g,S*ZS**2)+max(bmin*bmin,bmax*bmax)*Vtail;J2e=rows(cf(ge,S*ZS**2))+max(bmin*bmin,bmax*bmax)*et
 hg=cf(obj(o['HG'])[sl,sl],S**3);hge=rows(cf(obj(o['HGE'])[sl,sl],S**3));cx=cf(obj(o['C'])[sl,sl],S**3);ce=rows(cf(obj(o['CE'])[sl,sl]+obj(o['CE'])[sl,sl].T,S**3));Y2=cf(obj(o['Y2'])[sl,sl],S**4);Y2e=rows(cf(obj(o['Y2E'])[sl,sl],S**4));ht=rows(hg-Y2)+hge+Y2e;assert ht>=0;O=cx+cx.T-hg;oe=ce+hge+2*max(abs(bmin),bmax)*sqrtup(vt*ht);B2=J2-O;B2e=J2e+oe
 print(par,tag,'Berror',float(Be),'B2error',float(B2e),'vTail',float(vt),'exteriorError',float(oe),flush=True)
 return locals()

def audit(d,C,k,rot=None):
 # C maps from trial k to full eight columns. compare actual residual representation v0 before adding repr error.
 G=d['G'][:k,:k];Ge=rows(d['Ge'][:k,:k]);B=d['Bg'];B2=d['B2'];Be=d['Be'];B2e=d['B2e'];q=d['q'][:,:k];x=d['x'][:,:k];E=d['E'][:k,:k];diag=d['diag'][:k,:k];e=sqrtup(sum(v*v for v in d['er'][:k]));nq=sqrtup(rows(cf(q.T@q,S*S)));nx=norm(x,S*S);eE=2*nx*d['eK']*nq+nx*nx*d['eD'];derr=(d['eY']+sqrtup(d['tau']))*nq+(d['eU']+sqrtup(d['tauU']))*nx
 C=cf(C,S);nc=norm(np.vectorize(lambda v:int(v*S),otypes=[object])(C),S);BC=B[:k,:]@C;H=G-BC-BC.T+C.T@(B+B2)@C+diag;eh=Ge+2*nc*Be+nc*nc*(Be+B2e)
 # representation and operator-action errors in D^-1 norm; use <= I and combine all by triangle. Diag penalty included; this is a safe bound.
 h=rows(H)+eh;assert h>=0;dd=derr+e;allow=eE+(eh+2*sqrtup(h)*dd+dd*dd)/d['mu'];L=E-H/d['mu'];low=exact_lower(L)-allow;print(d['par'],d['tag'],'k',k,'head',float(h),'allow',float(allow),'lower',float(low),flush=True)
 if rot is not None:
  R=cf(rot,S);dr=rows(R.T@R-I(R.shape[1]));rl=exact_lower(R.T@L@R)-(1+dr)*allow;print('rot',R.shape[1],float(rl),flush=True)
 return L,allow,low
if __name__=='__main__':
 for par,tag in [('odd','surplus'),('even','raw'),('even','surplus')]:
  d=load(par,tag);B=d['Bg'];B2=d['B2'];M=np.array(B+B2,float);print('response eigenvalues',np.linalg.eigvalsh(M),flush=True)
  for k in [1,2,4,8]:
   C=np.rint(np.linalg.solve(M,np.array(B[:,:k],float))*S).astype(np.int64);L,allow,low=audit(d,C,k)
   if k==8:
    ev,R=np.linalg.eigh(np.array(L,float));print('central',ev,flush=True);R=np.rint(R*S).astype(np.int64)
    for n in range(1,9):
     r=R[:,-n:];rr=cf(r,S);dr=rows(rr.T@rr-I(n));rl=exact_lower(rr.T@L@rr)-(1+dr)*allow;print('rotated',n,float(rl),flush=True)
    np.savez_compressed(root/f'{par}_{tag}_response.npz',C=C,R=R,scale=S,lower=str(low),L=np.array(L,str),allow=str(allow))

def loadjoint(tag='joint'):
 d=load('odd','surplus');a=np.load(root/f'odd_{tag}_circle.npz');assert np.array_equal(a['ids'],d['c']['ids']);assert np.array_equal(obj(a['alpha'])[:len(d['alpha'])],d['alpha']);d['alpha']=obj(a['alpha']);d['Z']=obj(a['Z']);d['Ze']=obj(a['Ze']);d['G']=cf(obj(a['G']),S**3);d['Ge']=cf(obj(a['Ge']),S**3);d['er']=[F(str(v)) for v in a['er']];d['e']=sqrtup(sum(v*v for v in d['er']));d['t']=obj(a['t']);d['diag']=cf(d['t'].T@((10**8-d['alpha'])[:,None]*d['t']),10**8*S*S);zz,zze=gram(d['Z'],d['Ze']);d['Z2']=cf(zz,S**4);d['Z2E']=cf(zze,S**4);d['Vtail']=d['G']-d['Z2'];d['et']=rows(d['Ge']+d['Z2E']);d['vt']=rows(d['Vtail'])+d['et'];assert d['vt']>=0;return d
