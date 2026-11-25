from pathlib import Path
import sys, hashlib, json
from fractions import Fraction as F
workspace=Path(__file__).resolve().parent
par=sys.argv[1]; assert par in ('even','odd')
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
sys.path.insert(0,str(source))
from bv_joint_trial import *
code=(workspace/'cr_int_t1d_cb'/'cb_matrix.py').read_text();exec(code[:code.index("if __name__=='__main__':")])
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
assert par=='even' and '--fixed' not in sys.argv
assert np.load(source/'bu_trial.npz')['Q'].shape[1]==4
d=loadjoint();a=np.load(source/f'{par}_joint_Uaction.npz');b=build(d,a);n=len(b['E']);eye=b['eye'];Z=d['Z'];Ze=d['Ze'];U=obj(a['U']);Ue=obj(a['Ue']);Y=obj(a['Y']);Ye=obj(a['YE'])
m=b['mid'];r=b['rad'];s=m+F(1,2);c0=m*m+m+r*r
assert c0>=0 and r>=0
# Actual full u*h and v*h vanish. The finite Fourier synthesis centers below
# are rational, so divide exactly rather than casting fractional entries to int.
ww,wwe=gram(Z+Y-s*U,Ze+Ye+abs(s)*Ue)
Gu=b['Gu'];Gvu=b['Gvu'];full=d['G']+b['HG']+s*s*Gu-s*(Gvu+Gvu.T)
err=rows(d['Ge'])+b['HGe']+s*s*b['Gue']+2*abs(s)*b['Gvue']+rows(wwe/S**4)
W=full-ww/S**4+err*eye
assert np.array_equal(W,W.T)
c,ce=cross(U,Ue,Y,Ye,np.ones(len(U),object));midc=cf(c,S**4);mide=norm(ce,S**4)
base=b['base']-(d['bmax']+b['b2max'])*b['U']+c0*b['U']+m*(midc+midc.T)+2*abs(m)*mide*eye
Q=cf(d['q'].T@d['q'],S*S);X=cf(d['x'].T@d['x'],S**4);aq=d['eY']+sqrtup(d['tau']);ax=d['eU']+sqrtup(d['tauU']);er=d['e'];nq=sqrtup(rows(Q));nx=sqrtup(rows(X));v=[aq*nq,ax*nx,er];total=sum(v);lamb=[vi/total for vi in v];assert all(li>0 for li in lamb) and sum(lamb)==1
D2=aq**2/lamb[0]*Q+ax**2/lamb[1]*X+er**2/lamb[2]*eye
eta=nx/nq;assert eta>0
Edir=d['E']-d['eK']*(eta*Q+X/eta)-d['eD']*X
ql=rows(Q-eye)


from scipy.optimize import minimize
out=Path(sys.argv[2]);frozen=json.loads(Path(sys.argv[3]).read_text())
tail=json.loads((out/'contracted-tail.json').read_text())
assert tail['cutoff']==16384 and n==4
th,ga=map(F,frozen['parameters']);et=F(frozen['eta']);ws=list(map(F,frozen['weights']))
def compute(p,matrices):
 th,ga,et,ws=p
 assert min(th,ga,et)>0 and min(ws)>0 and sum(ws)==1 and len(ws)==len(matrices)
 H=base+r*(th*W+b['U']/th)
 D2=sum(A/w for A,w in zip(matrices,ws))
 L=d['E']-d['eK']*(et*Q+X/et)-d['eD']*X-((1+ga)*H+(1+1/ga)*D2)/d['mu']
 assert np.array_equal(L,L.T)
 return L
def publish(tag,p,mats):
 L=compute(p,mats);raw=exact_lower(L);low=raw/(1+ql) if raw>0 else raw
 rec={'n':n,'lower':str(low),'positive':bool(low>0),'parameters':[str(v) for v in p[:3]],'weights':[str(v) for v in p[3]],'matrix':[[str(v) for v in row] for row in L],'qdefect':str(ql)}
 if low>0:
  margin=F(int(low*100000//1),100000)
  if margin==low:margin-=F(1,100000)
  if margin>0:
   A=L-margin*(1+ql)*eye;T=np.vectorize(lambda v:int(v*S//1),otypes=[object])(A)-n*eye
   assert ldl_positive(cf(T,S));rec.update(margin=str(margin),upper_triangle=[int(T[i,j]) for i in range(n) for j in range(i,n)])
 (out/f'{tag}.json').write_text(json.dumps(rec,indent=2)+'\n')
 print('CJ EXACT',tag,rec['lower'],float(low),flush=True)
 return rec
baseline=publish('baseline',(th,ga,et,ws),[aq**2*Q,ax**2*X,er**2*eye])
assert baseline['matrix']==frozen['matrix'] and baseline['lower']==frozen['lower']
print('CH EXACT MATRIX MATCH',flush=True)
# Combined physical Neumann tail, after contraction, is bounded independently
# of the finite head-action and physical representation errors. Use triangle
# and weighted Cauchy-Schwarz; no assumption of inverse-metric orthogonality.
tt=list(map(F,tail['families']['combined']['norms']));assert len(tt)==n and all(v>=0 for v in tt)
tf=F(tail['families']['combined']['frobenius'])
# Frobenius and directional per-column alternatives.
families={'combined-uniform':[d['eY']**2*Q,d['eU']**2*X,tf**2*eye,er**2*eye]}
diag=[]
for i,t in enumerate(tt):
 A=np.zeros((n,n),object);A[i,i]=t*t;diag.append(A)
families['combined-directional']=[d['eY']**2*Q,d['eU']**2*X,*diag,er**2*eye]
# Retain old full-space certified bounds as a control; do not replace any
# inherited source/representation error except the tail.
Ef=np.asarray(d['E'],float);Bf=np.asarray(base,float);Wf=np.asarray(W,float);Uf=np.asarray(b['U'],float);Qf=np.asarray(Q,float);Xf=np.asarray(X,float);rf=float(r);mu=float(d['mu']);ek=float(d['eK']);ed=float(d['eD'])
def decode(z):
 p=np.exp(z[:3]);v=np.r_[z[3:],0.];v=np.exp(v-np.max(v));return p,v/v.sum()
res={'baseline':{k:v for k,v in baseline.items() if k!='matrix'}}
for tag,mats in families.items():
 a=[np.asarray(A,float) for A in mats]
 sizes=np.array([np.sqrt(max(1e-60,np.linalg.norm(A,2))) for A in a]);w=sizes/sizes.sum()
 z=np.r_[np.log([float(th),float(ga),float(et)]),np.log(w[:-1]/w[-1])]
 bounds=[(-9,5),(-9,5),(-12,12)]+[(-30,30)]*(len(w)-1)
 z=np.clip(z,[b[0] for b in bounds],[b[1] for b in bounds])
 def objective(y):
  (t,g,e),w=decode(y);H=Bf+rf*(t*Wf+Uf/t);D=sum(A/v for A,v in zip(a,w))
  return -float(np.linalg.eigvalsh(Ef-ek*(e*Qf+Xf/e)-ed*Xf-((1+g)*H+(1+1/g)*D)/mu)[0])
 opt=minimize(objective,z,method='Powell',bounds=bounds,options={'maxiter':150,'xtol':1e-8,'ftol':1e-12})
 p,w=decode(opt.x if objective(opt.x)<=objective(z) else z)
 aa=[max(1,int(v*10**12)) for v in w];weights=[F(v,sum(aa)) for v in aa]
 par=tuple(F(str(float(v))).limit_denominator(10**8) for v in p)+(weights,)
 rec=publish(tag,par,mats);res[tag]={k:v for k,v in rec.items() if k!='matrix'}
(out/'summary.json').write_text(json.dumps(res,indent=2)+'\n')
