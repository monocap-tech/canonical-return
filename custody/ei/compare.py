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
d=loadjoint()
# Loadjoint reads the historical 16384 metadata. Never use it for the extended
# finite representative. Authenticate and explicitly attach the new budgets.
root=workspace.parent.parent
eg=root/'custody/eg/generated';eh=root/'custody/eh/generated';dv=root/'custody/dv/generated';dz=root/'custody/dz/generated/trial31'
meta=json.loads((eh/'representation-budget.json').read_text())
ext=json.loads((eg/'extension.json').read_text())
assert meta['cutoff']==ext['cutoff']==65536
assert ext['head_half_target_met'] and ext['all_block_half_targets_met']
for k,j in [('head_error','head_error'),('head_tail','head_tail'),('physical_error','physical_finite_error'),('physical_tail','physical_tail'),('physical_blocks','physical_blocks')]:
    assert meta[k]==ext[j],k
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name in ('bs_trial.npz','bs_circle.npz'):
    assert digest(source/name)==digest(eh/name),name
assert digest(source/'odd_joint_Uaction.npz')==digest(dz/'odd_joint_Uaction.npz')
oldtrial=np.load(dz/'bs_trial.npz');newtrial=np.load(source/'bs_trial.npz')
assert np.array_equal(newtrial['ms'],np.arange(1,65536,2))
assert newtrial['yn'].shape==(32768,31)
assert np.array_equal(newtrial['yn'][:8192],oldtrial['yn'])
for k in oldtrial.files:
    if k not in ('yn','ms','errors'):assert np.array_equal(newtrial[k],oldtrial[k]),k
assert np.array_equal(d['c']['ids'],np.arange(1,32768,2))
head=np.load(eg/'odd_head_full.npz')
assert head['Y'].shape==(32768,119) and np.array_equal(head['indices'],newtrial['ms'])
assert F(str(head['error']))==F(meta['head_error'])==sum(map(F,map(str,head['tile_errors'])),F(0))
ht=np.load(workspace/'cr_int_t1d_bo'/'odd_tail.npz')
pt=np.load(dv/'merged/odd_merged_tail.npz')
assert F(str(ht['tail'][list(map(int,ht['cutoff'])).index(65536)]))==F(meta['head_tail'])
assert F(str(pt['tail'][list(map(int,pt['cutoff'])).index(65536)]))==F(meta['physical_tail'])
oldproj=np.load(dv/'merged/odd_combined_projection.npz')
assert F(str(oldproj['head_error']))==d['eY'] and F(str(oldproj['U_error']))==d['eU']
# Independently reconcile the exact new physical error from the disjoint chunks.
covered=[];newerr=F(0);common=F(0)
for chunk in ext['new_chunks']:
    lo,hi=chunk['lo'],chunk['hi'];p=np.load(eg/f'odd_U_{lo}_{hi}.npz')
    ids=p['indices'];assert np.array_equal(ids,np.arange(lo,hi+1,2)) and p['Y'].shape==(len(ids),352)
    assert int(p['scale'])==S
    ee=F(str(p['error']));assert ee==sum(map(F,map(str,p['tile_errors'])),F(0))
    extra=ee-F(len(ids)*352,S);assert extra>=0
    common+=extra;newerr+=ee;covered.extend(map(int,ids))
assert covered==list(range(16385,65536,2))
assert F(meta['physical_error'])==d['eU']+newerr
old_eU=d['eU'];old_eY=d['eY'];old_tau=d['tau'];old_tauU=d['tauU']
d.update(eY=F(meta['head_error']),eU=F(meta['physical_error']),tau=F(meta['head_tail']),tauU=F(meta['physical_tail']))
assert min(d[k] for k in ('eY','eU','tau','tauU'))>=0
d['derr']=(d['eY']+sqrtup(d['tau']))*d['nq']+(d['eU']+sqrtup(d['tauU']))*d['nx']
a=np.load(source/f'{par}_joint_Uaction.npz');b=build(d,a);n=len(b['E']);eye=b['eye'];Z=d['Z'];Ze=d['Ze'];U=obj(a['U']);Ue=obj(a['Ue']);Y=obj(a['Y']);Ye=obj(a['YE'])
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
assert frac(final_projection['U_error'])==old_eU
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
assert sum(b[3] for b in blocks)==old_eU and sum(b[4] for b in blocks)==d['tauU']
assert int(cutoffs[k])==65536
newblocks=[]
for (stage,lo,hi,e,t),v in zip(blocks,meta['physical_blocks']):
    assert (stage,lo,hi)==(v['stage'],v['lo'],v['hi'])
    enew=e+common+F(len(covered)*(hi-lo),S)
    assert F(v['finite_error'])==enew and F(v['tail_squared'])==t
    assert v['half_target_met']
    newblocks.append((stage,lo,hi,enew,t))
assert len(newblocks)==len(blocks)==19
blocks=newblocks
# Each block retains the entire common nonrounding reserve. Their error sum
# can exceed the aggregate bound; this is intentional, not an equality premise.
assert sum(b[3] for b in blocks)>=d['eU']
blockX=[cf(d['x'][lo:hi].T@d['x'][lo:hi],S**4) for _,lo,hi,_,_ in blocks]
assert np.array_equal(sum(blockX),X)
blocka=[e+sqrtup(t) for _,_,_,e,t in blocks]
assert all(v>0 for v in blocka)
print('EI EXACT PHYSICAL BLOCK DECOMPOSITION PASS',len(blocks),flush=True)
(out/'block-budget.json').write_text(json.dumps(dict(cutoff=str(cutoffs[k]),blocks=[dict(stage=s,lo=lo,hi=hi,error=str(e),tail=str(t),allowance=str(v)) for (s,lo,hi,e,t),v in zip(blocks,blocka)],scalar_allowance=str(ax)),indent=2)+'\n')
# Physical B is nonnegative, but the unrestricted circle symbol is signed.
# Its global lower endpoint is m-r (CB); the omitted U Gram is bounded by b['U'].
Bu_lower=b['Bu']-b['Bue']*eye+min(F(0),m-r)*b['U']
Cross_upper=b['Gvu']+b['Gvu'].T+2*b['Gvue']*eye
assert np.array_equal(Bu_lower,Bu_lower.T)
assert np.array_equal(Cross_upper,Cross_upper.T)
center=F(0)
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
    assert 0<=center<=F(1,2)
    # Let T=M^-1, e=V-MU.  T is self-adjoint with 0<=T<=I.
    # U+Te=(U+c e)+(T-cI)e; ||T-cI||<=1-c for 0<=c<=1/2.
    # Gram(U+c e)=(1-2c)Gram(U)+c sym(U*V)-2c U*BU+c^2 Gram(e).
    Ac=(1-2*center)*Gu_upper+center*Cross_upper-2*center*Bu_lower+center**2*R
    P=H-diag if mode==0 else (1+t)*Ac+(1+1/t)*(1-center)**2*R
    return E-(H+delta*P+(1+1/delta)*D)/d['mu'],H,R,P
def cj_publish(tag,par,mode):
    th,delta,t,et,ws=par
    assert all(v>0 for v in [th,delta,t,et,*ws]) and sum(ws)==1
    L,H,R,P=cj_matrix(*par,mode)
    assert all(np.array_equal(x,x.T) for x in [L,H,R,P])
    raw=exact_lower(L);low=raw/(1+ql) if raw>0 else raw
    rec=dict(mode=mode,center=str(center),n=n,split=len(ws)==len(blocks)+2,parameters=[str(v) for v in par[:4]],weights=[str(v) for v in ws],lower=str(low),positive=bool(low>0),qdefect=str(ql))
    for name,A in [('matrix',L),('H',H),('R',R),('P',P),('F_lower',F_lower),('Gu_upper',Gu_upper),('Bu_upper',Bu_upper),('Bu_lower',Bu_lower),('Cross_upper',Cross_upper)]:
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
    print('EI EXACT',tag,rec['lower'],float(low),rec['parameters'],rec['positive'],flush=True)
    return rec
# numerical proposals only, no float-derived bounds.
flF=np.asarray(F_lower,float);flG=np.asarray(Gu_upper,float);flDiag=np.asarray(diag,float)
blockaf=list(map(float,blocka));blockXf=[np.asarray(G,float) for G in blockX]
def cj_score(z,mode):
    p=np.exp(z[:4]);ww=np.exp(np.r_[z[4:],0]-max(0,*z[4:]));ww/=ww.sum()
    th,de,t,et=p
    H=Bf+rf*(th*Wf+Uf/th)
    R=H-flDiag-flF
    c=float(center)
    Ac=(1-2*c)*flG+c*np.asarray(Cross_upper,float)-2*c*np.asarray(Bu_lower,float)+c*c*R
    P=H-flDiag if mode==0 else (1+t)*Ac+(1+1/t)*(1-c)**2*R
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
# Evaluate EB's exact parameters first, then tune only the newly attached family.
oldbest=json.loads((root/'custody/eb/generated/result.json').read_text())['best']
assert oldbest['lower']=='-3/5000' and oldbest['center']=='1/2'
start=(*map(F,oldbest['parameters']),list(map(F,oldbest['weights'])))
assert len(start[-1])==21
recs=[]
for c in (F(0),F(1,4),F(1,2)):
    center=c
    recs.append(cj_publish('matched-eb-parameters-center-'+str(c).replace('/','_'),start,1))
# Rebalance the weights to the actual reduced representation allowances.
a0=[aq*nq]+[v*sqrtup(rows(G)) for v,G in zip(blocka,blockX)]+[er]
assert all(v>0 for v in a0)
s0=sum(a0);newws=[v/s0 for v in a0]
starts=[start,(*start[:4],newws)]
trials=[]
bounds=[(-9,5),(-12,5),(-9,9),(-12,12)]+[(-40,40)]*20
for c in (F(0),F(1,4),F(1,2)):
    center=c
    for i,p in enumerate(starts):
        z=np.r_[np.log([float(v) for v in p[:4]]),np.log([float(v/p[-1][-1]) for v in p[-1][:-1]])]
        z=np.clip(z,[v[0] for v in bounds],[v[1] for v in bounds])
        opt=minimize(lambda x:-cj_score(x,1),z,method='L-BFGS-B',bounds=bounds,options={'maxiter':150})
        if np.all(np.isfinite(opt.x)) and np.isfinite(opt.fun) and cj_score(opt.x,1)>=cj_score(z,1):z=opt.x
        recs.append(cj_publish('optimized-'+str(c).replace('/','_')+'-'+str(i),cj_rational(z),1))
        trials.append((cj_score(z,1),c,z))
_,center,z=max(trials,key=lambda x:x[0])
opt=minimize(lambda x:-cj_score(x,1),z,method='Powell',bounds=bounds,options={'maxiter':150})
if np.all(np.isfinite(opt.x)) and np.isfinite(opt.fun) and cj_score(opt.x,1)>=cj_score(z,1):
    recs.append(cj_publish('optimized-Powell',cj_rational(opt.x),1))
best=max(recs,key=lambda x:F(x['lower']))
omit=['matrix','H','R','P','F_lower','Gu_upper','Bu_upper','Bu_lower','Cross_upper']
summary=[{k:v for k,v in rec.items() if k not in omit} for rec in recs]
report=dict(scope='extended physical cutoff, matched circular cutoff; inherited analytical enclosures',dimension=n,physical_correctors=352,physical_cutoff=65536,circular_cutoff=32768,head_allowance=str(aq),physical_allowance=str(ax),representation_error=str(er),old_lower=oldbest['lower'],comparisons=summary,best={k:v for k,v in best.items() if k not in omit})
(out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print('EI BEST',best['center'],best['lower'],best['positive'],flush=True)
