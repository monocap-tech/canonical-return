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

# CI is a diagnostic, not an acceptance promotion. Reproduce CH exactly first.
from math import isqrt
assert not frozen['split']
th,ga=map(F,frozen['parameters']);et=F(frozen['eta']);ws=list(map(F,frozen['weights']))
p=(th,ga,et,ws)
rec=evaluate(p)
assert rec['matrix']==frozen['matrix'] and rec['lower']==frozen['lower'] and rec['qdefect']==frozen['qdefect']
print('CH EXACT MATRIX MATCH',rec['lower'],flush=True)
M=np.asarray([[F(x) for x in row] for row in rec['matrix']],dtype=object)
mu=d['mu'];assert mu>0
fac=(1+ga)/mu;dfac=(1+1/ga)/mu
pieces={
 'mixed_energy_Q':d['eK']*et*Q,
 'mixed_energy_X':d['eK']/et*X,
 'energy_D':d['eD']*X,
 'residual_base':fac*base,
 'omitted_shifted_gram':fac*r*th*W,
 'omitted_inverse_gram':fac*r*b['U']/th,
 'source_Y_and_tail':dfac*aq**2/ws[0]*Q,
 'source_U_and_tail':dfac*ax**2/ws[1]*X,
 'representation':dfac*er**2/ws[2]*eye,
}
assert np.array_equal(d['E']-sum(pieces.values()),M)
# Float eigendirection is merely a proposal. All contractions below exact.
eig,vec=np.linalg.eigh(np.asarray(M,float))
z=np.array([F(int(x)) for x in np.rint(vec[:,0]*10**8)],dtype=object)
zz=z@z;assert zz>0
def ray(A):
 return (z@A@z)/zz
def sqlo(x):
 x=F(x);assert x>=0
 s=10**40
 t=F(isqrt((x.numerator*s*s)//x.denominator),s)
 assert t>=0 and t*t<=x
 return t
vals={k:ray(v) for k,v in pieces.items()}
assert ray(d['E'])-sum(vals.values())==ray(M)
# Remove individual pieces only as counterfactual diagnostics, never as proof.
leave={}
for k,A in pieces.items():
 leave[k]={'ray':str(ray(M+A)), 'float_min_eig':float(np.linalg.eigvalsh(np.asarray(M+A,float))[0])}
# Uniform or per-column error allocation. For ANY positive weights summing to
# one, Cauchy-Schwarz gives sum a_i^2/lambda_i >= (sum a_i)^2.
# For ANY positive theta, theta*w+u/theta >= 2 sqrt(w*u).
# Similarly for eta and gamma. All square roots rounded DOWN for a rigorous
# upper bound on the comparator along z. No numerical optimizer needed.
q,x,w,u=map(ray,[Q,X,W,b['U']])
assert min(q,x,w,u)>=0 and min(aq,ax,er,r,d['eK'],d['eD'])>=0
hlo=ray(base)+2*r*sqlo(w*u)
assert hlo>=0, 'nonnegative residual lower estimate required'
ylo=aq*sqlo(q);ulo=ax*sqlo(x)
dlo=(ylo+ulo+er)**2
splitlo=(ylo+ulo+sum(ei*abs(zi) for ei,zi in zip(d['er'],z))/ (sqlo(zz)+F(1,10**40)))**2
# sqlo(zz)+grid is an upper norm bound; use it to lower-bound split term.
assert all(ei>=0 for ei in d['er'])
elo=ray(d['E'])-2*d['eK']*sqlo(q*x)-d['eD']*x
def universal(dl):
 return elo-(hlo+dl+2*sqlo(hlo*dl))/mu
bounds={'uniform':universal(dlo),'split':universal(splitlo)}
# These bounds only constrain this family of sufficient lower matrices.
# Negative value does not give a negative direction of the original form.
report={
 'baseline':rec,'direction':[str(v) for v in z],'direction_norm_sq':str(zz),
 'rayleigh':{k:str(v) for k,v in vals.items()},
 'rayleigh_decimal':{k:float(v) for k,v in vals.items()},
 'center_energy':str(ray(d['E'])),'center_energy_decimal':float(ray(d['E'])),
 'comparison_ray':str(ray(M)),'comparison_ray_decimal':float(ray(M)),
 'leave_one_out_diagnostics':leave,
 'universal_upper_ray':{k:str(v) for k,v in bounds.items()},
 'universal_upper_ray_decimal':{k:float(v) for k,v in bounds.items()},
 'scalar_family_excluded':{k:bool(v<0) for k,v in bounds.items()},
 'coercivity':str(mu),'H_direction_lower':str(hlo),
 'D_direction_lower_uniform':str(dlo),'D_direction_lower_split':str(splitlo),
 'energy_direction_upper':str(elo),
 'matrices':{k:[[str(v) for v in row] for row in A] for k,A in pieces.items()}
}
(out/'diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
print('CI RAYLEIGH',json.dumps({'center':float(ray(d['E'])),'comparison':float(ray(M)),'terms':report['rayleigh_decimal'],'mu':float(mu)}),flush=True)
print('CI UNIVERSAL',json.dumps({'upper':report['universal_upper_ray_decimal'],'excluded':report['scalar_family_excluded']}),flush=True)
# Preserve the actual source modules used for later structural refinement.
import shutil
for pth in [source/'bv_joint_trial.py',workspace/'cr_int_t1d_cb'/'cb_matrix.py',source/'bv_joint_action.py']:
 shutil.copyfile(pth,out/pth.name)
