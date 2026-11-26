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
assert par=='odd' and '--fixed' not in sys.argv
shape=np.load(source/'bs_trial.npz')['Q'].shape
assert shape[0]==119 and shape[1]==31
assert np.load(source/'bs_trial.npz')['X'].shape==(352,119)
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
frozen=None
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
 print('CJ CONTROL',tag,rec['lower'],float(F(rec['lower'])),rec['parameters'],rec['positive'],flush=True)
 return rec

# CB historical scalar-tail comparison, not the later secant comparison.
if '--control' in sys.argv:
    old=np.load(source/'cb_matrix_accepted.npz')
    pp=list(map(F,old['parameters']))
    L,H=compare(b,*pp)
    raw=exact_lower(L);low=raw/(1+ql) if raw>0 else raw
    assert str(low)==str(old['lower'])
    assert np.array_equal(np.array(L,str),old['L']) and np.array_equal(np.array(H,str),old['H'])
    assert low>0
    print('CB EXACT MATRIX MATCH',str(low),flush=True)
    (out/'cb-control.json').write_text(json.dumps(dict(lower=str(low),parameters=list(map(str,pp)),matrix_match=True))+'\n')
bp=(F(1,50),F(1,175),eta,lamb)
rec=publish('secant-control',bp)
Ef=np.asarray(d['E'],float);Bf=np.asarray(base,float);Wf=np.asarray(W,float);Uf=np.asarray(b['U'],float);Qf=np.asarray(Q,float);Xf=np.asarray(X,float);If=np.eye(n);rf=float(r);mu=float(d['mu']);ek=float(d['eK']);ed=float(d['eD'])
def rational_weights(w):
    a=[max(1,int(x*10**10)) for x in w]
    return [F(x,sum(a)) for x in a]
# H bounds 2<V,U>-<U,MU>+||V-MU||^2 plus the fixed diagonal.
# F_lower bounds the variational term from below using actual physical Grams.
# B_u is the enclosed finite symbol-weighted Gram plus the omitted tail.
Gu_upper=b['Gu']+b['Gue']*eye
Bu_upper=b['Bu']+b['Bue']*eye+d['bmax']*b['U']
F_lower=b['Gvu']+b['Gvu'].T-2*b['Gvue']*eye-Gu_upper-Bu_upper
assert np.array_equal(F_lower,F_lower.T)
diag=d['diag']
assert d['mu']>0
# Recover the disjoint physical blocks from the original merge custody.
root=workspace.parent.parent
stages=('cm','co','cq','cr','ct','cv','cy','da','db','dc','de','dg','di','dk','dn','ds','du','dv')
def load(stage,rel):
    return np.load(root/'custody'/stage/'generated'/rel)
def frac(x):return F(str(x))
assert d['x'].shape==(352,n)
blocks=[]
last_e=None;last_tau=None
final_projection=load('dv','merged/odd_combined_projection.npz')
final_tail=load('dv','merged/odd_merged_tail.npz')
cutoffs=final_tail['cutoff']
selected=[i for i,v in enumerate(final_tail['tail']) if frac(v)==d['tauU']]
assert len(selected)==1,('tail identity',list(map(str,cutoffs)),str(d['tauU']))
k=selected[0]
assert frac(final_projection['U_error'])==d['eU']
for stage in stages:
    u=load(stage,'new/odd_U_full.npz')
    v=load(stage,'new/odd_U_tail.npz')
    p=load(stage,'merged/odd_combined_projection.npz')
    t=load(stage,'merged/odd_merged_tail.npz')
    e=frac(u['error']);tau=list(map(frac,v['tail']))
    total_e=frac(p['U_error']);total_tau=list(map(frac,t['tail']))
    assert np.array_equal(v['cutoff'],cutoffs) and np.array_equal(t['cutoff'],cutoffs)
    assert u['Y'].shape==(8192,16) and int(u['scale'])==int(p['Yscale'])==int(final_projection['Yscale'])
    assert e==sum(map(frac,u['tile_errors']),F(0))
    if last_e is None:
        last_e=total_e-e;last_tau=[x-y for x,y in zip(total_tau,tau)]
        blocks.append(('base64',0,64,last_e,last_tau[k]))
        assert np.array_equal(p['Y'][:,:183],final_projection['Y'][:,:183])
    assert last_e+e==total_e and [x+y for x,y in zip(last_tau,tau)]==total_tau
    lo=64+16*(len(blocks)-1);hi=lo+16
    assert p['Y'].shape==(8192,119+hi)
    assert np.array_equal(u['Y'],p['Y'][:,-16:])
    assert np.array_equal(u['Y'],final_projection['Y'][:,119+lo:119+hi])
    blocks.append((stage,lo,hi,e,tau[k]))
    last_e=total_e;last_tau=total_tau
assert blocks[-1][2]==352 and len(blocks)==19
assert all(e>=0 and t>=0 for _,_,_,e,t in blocks)
assert sum(b[3] for b in blocks)==d['eU'] and sum(b[4] for b in blocks)==d['tauU']
blockX=[cf(d['x'][lo:hi].T@d['x'][lo:hi],S**4) for _,lo,hi,_,_ in blocks]
assert np.array_equal(sum(blockX),X)
blocka=[e+sqrtup(t) for _,_,_,e,t in blocks]
assert all(v>0 for v in blocka)
print('DW EXACT PHYSICAL BLOCK DECOMPOSITION PASS',len(blocks),flush=True)
(out/'block-budget.json').write_text(json.dumps(dict(cutoff=str(cutoffs[k]),blocks=[dict(stage=s,lo=lo,hi=hi,error=str(e),tail=str(t),allowance=str(v)) for (s,lo,hi,e,t),v in zip(blocks,blocka)],scalar_allowance=str(ax)),indent=2)+'\n')
def cj_matrix(th,delta,t,et,ws,mode):
    H=base+r*(th*W+b['U']/th)
    R=H-diag-F_lower
    assert len(ws) in (3,len(blocks)+2)
    if len(ws)==3:
        D=aq**2/ws[0]*Q+ax**2/ws[1]*X+er**2/ws[2]*eye
    else:
        # Triangle and weighted Cauchy-Schwarz, one weight per physical block.
        D=aq**2/ws[0]*Q+er**2/ws[-1]*eye
        for v,G,w in zip(blocka,blockX,ws[1:-1]): D=D+v*v/w*G
    E=d['E']-d['eK']*(et*Q+X/et)-d['eD']*X
    # mode 0: charge representation interaction only against physical inverse energy.
    # mode 1: use M^-2 via U + M^-1(V-MU); M>=I inherited from BU.
    P=H-diag if mode==0 else (1+t)*Gu_upper+(1+1/t)*R
    return E-(H+delta*P+(1+1/delta)*D)/d['mu'],H,R,P
def cj_publish(tag,par,mode):
    th,delta,t,et,ws=par
    assert all(v>0 for v in [th,delta,t,et,*ws]) and sum(ws)==1
    L,H,R,P=cj_matrix(*par,mode)
    assert all(np.array_equal(x,x.T) for x in [L,H,R,P])
    raw=exact_lower(L);low=raw/(1+ql) if raw>0 else raw
    rec=dict(mode=mode,n=n,split=len(ws)==len(blocks)+2,parameters=[str(v) for v in par[:4]],weights=[str(v) for v in ws],lower=str(low),positive=bool(low>0),qdefect=str(ql))
    for name,A in [('matrix',L),('H',H),('R',R),('P',P),('F_lower',F_lower),('Gu_upper',Gu_upper),('Bu_upper',Bu_upper)]:
        rec[name]=[[str(v) for v in row] for row in A]
    if low>0:
        margin=F(int(low*100000//1),100000)
        if margin==low:margin-=F(1,100000)
        if margin>0:
            A=L-margin*(1+ql)*eye
            T=np.vectorize(lambda v:int(v*S//1),otypes=[object])(A)-n*eye
            assert ldl_positive(cf(T,S))
            rec.update(margin=str(margin),upper_triangle=[int(T[i,j]) for i in range(n) for j in range(i,n)])
    (out/f'{tag}.json').write_text(json.dumps(rec,indent=2)+'\n')
    print('DW EXACT',tag,rec['lower'],float(low),rec['parameters'],rec['positive'],flush=True)
    return rec
# numerical proposals only, no float-derived bounds.
flF=np.asarray(F_lower,float);flG=np.asarray(Gu_upper,float);flDiag=np.asarray(diag,float)
blockaf=list(map(float,blocka));blockXf=[np.asarray(G,float) for G in blockX]
def cj_score(z,mode):
    p=np.exp(z[:4]);ww=np.exp(np.r_[z[4:],0]-max(0,*z[4:]));ww/=ww.sum()
    th,de,t,et=p
    H=Bf+rf*(th*Wf+Uf/th)
    R=H-flDiag-flF
    P=H-flDiag if mode==0 else (1+t)*flG+(1+1/t)*R
    if len(ww)==3:
        D=float(aq)**2/ww[0]*Qf+float(ax)**2/ww[1]*Xf+float(er)**2/ww[2]*If
    else:
        D=float(aq)**2/ww[0]*Qf+float(er)**2/ww[-1]*If
        for v,G,w in zip(blockaf,blockXf,ww[1:-1]):D=D+v*v/w*G
    L=Ef-ek*(et*Qf+Xf/et)-ed*Xf-(H+de*P+(1+1/de)*D)/mu
    return float(np.linalg.eigvalsh(L)[0])
def cj_rational(z):
    p=np.exp(z[:4]);ww=np.exp(np.r_[z[4:],0]-max(0,*z[4:]));ww/=ww.sum()
    return (*[F(str(float(x))).limit_denominator(10**7) for x in p],rational_weights(ww))
# Reproduce frozen baseline before experimenting.
frozen_dir=Path(sys.argv[3])
old=json.loads((frozen_dir/'cj-mode1-start1.json').read_text())
p=(*map(F,old['parameters']),list(map(F,old['weights'])))
control=cj_publish('scalar-control',p,1)
for key in ('matrix','H','R','P','F_lower','Gu_upper','Bu_upper','lower','qdefect'):
    assert control[key]==old[key],key
print('DW EXACT FROZEN SCALAR CONTROL MATCH',flush=True)
# Diagnose the limiting eigenvector of the frozen exact comparison.
# These floating diagnostics are not certificate premises.
th,de,t,et,ws=p
L,H,R,P=cj_matrix(*p,1)
eig,Vlim=np.linalg.eigh(np.asarray(L,float));zlim=Vlim[:,0]
def ray(A):return float(zlim@np.asarray(A,float)@zlim)
E=d['E']-d['eK']*(et*Q+X/et)-d['eD']*X
pieces=dict(central_energy=ray(d['E']),energy_cross_reserve=ray(d['eK']*(et*Q+X/et)),energy_physical_reserve=ray(d['eD']*X),H=ray(H/d['mu']),inverse_response=ray(de*P/d['mu']),head_representation=ray((1+1/de)*aq**2/ws[0]*Q/d['mu']),physical_representation=ray((1+1/de)*ax**2/ws[1]*X/d['mu']),column_representation=ray((1+1/de)*er**2/ws[2]*eye/d['mu']))
recon=pieces['central_energy']-sum(v for k,v in pieces.items() if k!='central_energy')
assert abs(recon-float(eig[0]))<1e-8
diagnostic=dict(smallest_eigenvalue=float(eig[0]),rayleigh_terms=pieces,rayleigh_reconstruction=recon,aq=str(aq),ax=str(ax),er=str(er),eU=str(d['eU']),tauU=str(d['tauU']),mu=str(d['mu']),blocks=[dict(stage=b[0],gram_rayleigh=ray(G),allowance=float(v),weighted_norm=float(v)*np.sqrt(max(0,ray(G)))) for b,G,v in zip(blocks,blockX,blocka)],scope='floating diagnostic only')
(out/'limiting-direction.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
print('DW LIMITING DIRECTION',diagnostic,flush=True)
# Freeze scalar comparison; distribute its physical weight among blocks.
th,de,t,et,ws=p
alloc=[v*v for v in blocka]
den=sum(alloc)
splitws=[ws[0]]+[ws[1]*v/den for v in alloc]+[ws[2]]
start=(th,de,t,et,splitws)
assert sum(splitws)==1
recs=[control,cj_publish('split-initial',start,1)]
z=np.r_[np.log([float(v) for v in start[:4]]),np.log([float(v/splitws[-1]) for v in splitws[:-1]])]
bounds=[(-9,5),(-9,5),(-9,9),(-12,12)]+[(-40,40)]*(len(splitws)-1)
z=np.clip(z,[v[0] for v in bounds],[v[1] for v in bounds])
for method in ('L-BFGS-B','Powell'):
    opt=minimize(lambda x:-cj_score(x,1),z,method=method,bounds=bounds,options={'maxiter':150})
    if np.all(np.isfinite(opt.x)) and np.isfinite(opt.fun):
        recs.append(cj_publish('split-'+method,cj_rational(opt.x),1))
        if -opt.fun>=cj_score(z,1):z=opt.x
best=max(recs,key=lambda x:F(x['lower']))
summary=[{k:v for k,v in rec.items() if k not in ['matrix','H','R','P','F_lower','Gu_upper','Bu_upper']} for rec in recs]
(out/'result.json').write_text(json.dumps(dict(dimension=n,physical_correctors=352,scalar_control=summary[0],comparisons=summary,best={k:v for k,v in best.items() if k not in ['matrix','H','R','P','F_lower','Gu_upper','Bu_upper']}),indent=2)+'\n')
print('DW BEST',best['split'],best['lower'],flush=True)
