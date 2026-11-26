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
 assert par=='even'
 a=np.load(root/'bu_trial.npz');c=np.load(root/'joint_circle.npz');b=np.load(parent/'cr_int_t1d_bu'/'bu_symbol.npz');en=np.load(parent/'cr_int_t1d_bu'/'even_BQ_energy.npz');q=obj(a['Q']);X=obj(a['X']);x=X@q;t=obj(a['t']);alpha=obj(a['alpha']);eY=F(str(a['head_error']));eU=F(str(a['Uerror']));tau=F(str(a['tail']));tauU=F(str(a['Utail']));mu=F(str(a['mu']));eps=F(str(a['eps']));sl=slice(None)

 K=obj(en['K']);D=obj(en['D']);KS=int(en['scale']);eK=F(str(en['Kerror']));eD=F(str(en['Derror']));nq=sqrtup(rows(cf(q.T@q,S*S)));nx=norm(x,S*S);E=(1-eps)*cf(q.T@q,S*S)+cf((x.T@K@q+q.T@K.T@x)*S+x.T@D@x,KS*S**4);eE=2*nx*eK*nq+nx*nx*eD;derr=(eY+sqrtup(tau))*nq+(eU+sqrtup(tauU))*nx
 diag=cf(t.T@((10**8-alpha)[:,None]*t),10**8*S*S)
 Z=obj(c['Z'])[:,sl];Ze=obj(c['Zerr'])[:,sl];ZS=int(c['Zscale']);G=cf(obj(c['G'])[sl,sl],int(c['Gscale']));Ge=cf(obj(c['Gerr'])[sl,sl],int(c['Gscale']));er=[F(str(e)) for e in c['reprerr'][sl]];e=sqrtup(sum(v*v for v in er));assert np.array_equal(c['ids'],b['ids']);WS=int(b['scale']);bl=[F(6,5)*(F(399,50)-F(int(v),WS)) for v in b['au']];bh=[F(6,5)*(F(399,50)-F(int(v),WS)) for v in b['al']];bmin=F(str(b['bmin']));bmax=F(6,5)*(F(399,50)+F(str(b['sumupper'])));bn=F(2394,125);assert all(l<=h for l,h in zip(bl,bh));bmin=min(bmin,min(bl));bmax=max(bmax,max(bh));
 # actual B gram center, with global signed-symbol tail interval rather than assume sign. Tail is G-Z*Z.
 zz,zze=gram(Z,Ze);Z2=cf(zz,ZS**2);Z2E=cf(zze,ZS**2);Vtail=G-Z2;et=rows(Ge+Z2E);vt=rows(Vtail)+et;assert vt>=0
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

def loadjoint(tag="joint"):
 return load("even","joint")
