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

# CF's stored comparison is the exact baseline. All later proposals are checked
# using rational arithmetic; floats can only select positive weights.
from scipy.optimize import minimize
out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
frozen=json.loads(Path(sys.argv[3]).read_text())
def evaluate(p, split=False):
 th,ga,et,weights=p
 assert th>0 and ga>0 and et>0 and all(w>0 for w in weights) and sum(weights)==1
 if split:
  assert len(weights)==n+2
  dd=aq**2/weights[0]*Q+ax**2/weights[1]*X+np.diag([v*v/w for v,w in zip(d['er'],weights[2:])]).astype(object)
 else:
  assert len(weights)==3
  dd=aq**2/weights[0]*Q+ax**2/weights[1]*X+er**2/weights[2]*eye
 hh=base+r*(th*W+b['U']/th)
 ee=d['E']-d['eK']*(et*Q+X/et)-d['eD']*X
 M=ee-((1+ga)*hh+(1+1/ga)*dd)/d['mu']
 assert np.array_equal(M,M.T)
 raw=exact_lower(M);low=raw/(1+ql) if raw>0 else raw
 rec=dict(n=n,parity=par,parameters=[str(th),str(ga)],eta=str(et),weights=[str(w) for w in weights],split=split,lower=str(low),positive=bool(low>0),matrix=[[str(x) for x in row] for row in M],qdefect=str(ql))
 if low>0:
  margin=F(int(low*100000//1),100000)
  if margin==low: margin-=F(1,100000)
  if margin>0:
   A=M-margin*(1+ql)*eye
   T=np.vectorize(lambda v:int(v*S//1),otypes=[object])(A)-n*eye
   assert ldl_positive(cf(T,S))
   rec.update(margin=str(margin),upper_triangle=[int(T[i,j]) for i in range(n) for j in range(i,n)])
 return rec
def publish(tag,p,split=False):
 rec=evaluate(p,split)
 (out/f'{tag}.json').write_text(json.dumps(rec,indent=2)+'\n')
 print('CI EXACT',tag,rec['lower'],float(F(rec['lower'])),rec['parameters'],rec['positive'],flush=True)
 return rec
assert not frozen['split']
bp=(*map(F,frozen['parameters']),F(frozen['eta']),list(map(F,frozen['weights'])))
rec=publish('cf-baseline',bp)
if '--control' in sys.argv:
 assert rec['lower']==frozen['lower'] and rec['matrix']==frozen['matrix'] and rec['qdefect']==frozen['qdefect']
 print('CH EXACT MATRIX MATCH',flush=True)
# Independent exact evaluation after parameter selection. Bound search is
# heuristic only, not an exhaustion or global optimum claim.
Ef=np.asarray(d['E'],float);Bf=np.asarray(base,float);Wf=np.asarray(W,float);Uf=np.asarray(b['U'],float);Qf=np.asarray(Q,float);Xf=np.asarray(X,float);If=np.eye(n);rf=float(r);mu=float(d['mu']);ek=float(d['eK']);ed=float(d['eD'])
def fl(th,ga,et,ws,split):
 if split: D=float(aq)**2/ws[0]*Qf+float(ax)**2/ws[1]*Xf+np.diag([float(v)**2/w for v,w in zip(d['er'],ws[2:])])
 else: D=float(aq)**2/ws[0]*Qf+float(ax)**2/ws[1]*Xf+float(er)**2/ws[2]*If
 H=Bf+rf*(th*Wf+Uf/th)
 return Ef-ek*(et*Qf+Xf/et)-ed*Xf-((1+ga)*H+(1+1/ga)*D)/mu
def score(th,ga,et,ws,split=False):
 return float(np.linalg.eigvalsh(fl(th,ga,et,ws,split))[0])
best=(score(*[float(v) for v in bp[:3]],list(map(float,bp[3]))),bp)
for ti in range(1,151):
 th=F(1,ti)
 for ga in sorted(set([F(1,j) for j in range(1,321)]+[F(j,10) for j in range(1,31)])):
  ev=score(float(th),float(ga),float(eta),list(map(float,lamb)))
  if ev>best[0]:best=(ev,(th,ga,eta,lamb))
expanded=publish('expanded-grid',best[1])
# Every weight is positive and exactly sums to one. Split case uses the
# per-column representation radii rather than their uniform Euclidean bound.
def rational_weights(w):
 scale=10**10
 a=[max(1,int(x*scale)) for x in w]
 return [F(x,sum(a)) for x in a]
def decode(z):
 p=np.exp(z[:3]);v=np.concatenate([z[3:],np.zeros(1)]);v=np.exp(v-np.max(v));return p,v/v.sum()
def rationalize(z):
 p,w=decode(z)
 return tuple(F(str(float(x))).limit_denominator(10**6) for x in p)+(rational_weights(w),)
def refine(initial,split):
 p,w=initial[:3],np.asarray(initial[3],float)
 z=np.concatenate([np.log([float(x) for x in p]),np.log(w[:-1]/w[-1])])
 bounds=[(-9,5),(-9,5),(-12,12)]+[(-25,25)]*(len(w)-1)
 z=np.clip(z,[b[0] for b in bounds],[b[1] for b in bounds])
 def objective(y):
  v,l=decode(y)
  return -score(*v,l,split)
 opt=minimize(objective,z,method='Powell',bounds=bounds,options={'maxiter':150,'xtol':1e-8,'ftol':1e-12})
 print('optimizer diagnostic',split,opt.success,opt.nfev,opt.fun,flush=True)
 candidate=rationalize(opt.x)
 # Never discard the starting exact comparison simply because a float search moved.
 if objective(opt.x)>objective(z): candidate=rationalize(z)
 return candidate
p=refine(best[1],False);joint=publish('joint-weights',p)
# Triangle inequality gives representation norm <= sum_i er_i |z_i|;
# weighted Cauchy-Schwarz gives the exact diagonal majorant in evaluate().
er2=[v*v for v in d['er']];assert all(v>=0 for v in er2)
sw=list(p[3][:2])+[p[3][2]*(v+F(1,10**80))/(sum(er2)+n*F(1,10**80)) for v in er2]
assert sum(sw)==1
ps=refine((*p[:3],sw),True);split=publish('split-representation',ps,True)
results={q['split'].__str__()+str(i):{k:v for k,v in q.items() if k!='matrix'} for i,q in enumerate([rec,expanded,joint,split])}
(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
print('CI SUMMARY',json.dumps(results),flush=True)

# Directional diagnostic only: float eigenvector is not a certified lower bound.
# Store exact component Rayleigh quotients in a rationalized test direction.
p=(*map(F,split['parameters']),F(split['eta']),list(map(F,split['weights'])))
th,ga,et,ws=p
M=fl(float(th),float(ga),float(et),list(map(float,ws)),True)
_,eigvec=np.linalg.eigh(M)
z=np.array([F(int(round(v*10**9)),10**9) for v in eigvec[:,0]],object)
assert z@z>0
pieces={
 'central_energy':d['E'],
 'energy_cross_error':-d['eK']*(et*Q+X/et),
 'energy_quadratic_error':-d['eD']*X,
 'residual_base':-(1+ga)*base/d['mu'],
 'omitted_shifted_v_tail':-(1+ga)*r*th*W/d['mu'],
 'omitted_u_tail':-(1+ga)*r*b['U']/th/d['mu'],
 'residual_q_error':-(1+1/ga)*aq**2/ws[0]*Q/d['mu'],
 'residual_x_error':-(1+1/ga)*ax**2/ws[1]*X/d['mu'],
 'residual_representation_error':-(1+1/ga)*np.diag([v*v/w for v,w in zip(d['er'],ws[2:])]).astype(object)/d['mu']}
actual=np.array([[F(v) for v in row] for row in split['matrix']],object)
assert np.array_equal(sum(pieces.values()),actual)
quot={k:str((z@v@z)/(z@z)) for k,v in pieces.items()}
assert sum(map(F,quot.values()))==(z@actual@z)/(z@z)
recdiag=dict(direction=[str(v) for v in z],terms=quot,rayleigh=str(sum(map(F,quot.values()))),diagnostic_only=True,action_rows=len(a['ids']),interior=len(a['j']),exterior=len(a['idx']))
(out/'directional-loss.json').write_text(json.dumps(recdiag,indent=2)+'\n')
print('CI DIRECTIONAL LOSS',json.dumps({k:float(F(v)) for k,v in quot.items()}), 'rayleigh',float(F(recdiag['rayleigh'])),flush=True)
