"""Independent receiver audit of completed EI evidence.
No import of the producer's comparison, synthesis, or LDL routines.
Requires an authenticated completed EI evidence checkout.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
import numpy as np
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parents[2];c=root/'custody';out=c/'ej/generated';out.mkdir(parents=True,exist_ok=True)
S=10**12;n=31
def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def auth(p):
    for line in (p/'sha256.txt').read_text().splitlines():
        h,name=line.split(maxsplit=1)
        assert digest(p/name.removeprefix('./'))==h,name
for stage in ['dv','dz','eb','eg','eh','ei']:auth(c/stage/'generated')
eg=c/'eg/generated';eh=c/'eh/generated';ei=c/'ei/generated';dv=c/'dv/generated';dz=c/'dz/generated/trial31'
def read(p):return json.loads(p.read_text())
ext=read(eg/'extension.json');meta=read(eh/'representation-budget.json');result=read(ei/'result.json')
assert ext['cutoff']==meta['cutoff']==result['physical_cutoff']==65536
assert result['circular_cutoff']==32768 and result['dimension']==n and result['physical_correctors']==352
assert ext['head_half_target_met'] and ext['all_block_half_targets_met']
for k,j in [('head_error','head_error'),('head_tail','head_tail'),('physical_error','physical_finite_error'),('physical_tail','physical_tail'),('physical_blocks','physical_blocks')]:assert meta[k]==ext[j]
old=np.load(dz/'bs_trial.npz');a=np.load(eh/'bs_trial.npz');circle=np.load(eh/'bs_circle.npz')
assert np.array_equal(a['ms'],np.arange(1,65536,2)) and a['yn'].shape==(32768,n)
assert np.array_equal(a['yn'][:8192],old['yn'])
for k in old.files:
    if k not in ['yn','ms','errors']:assert np.array_equal(a[k],old[k]),k
assert int(a['scale'])==S
q=a['Q'].astype(object);x=a['X'].astype(object)@q
assert q.shape==(119,n) and x.shape==(352,n)
h=np.load(eg/'odd_head_full.npz');pr=np.load(dv/'merged/odd_combined_projection.npz')
assert np.array_equal(h['Y'][:8192],pr['Y'][:,:119])
assert np.array_equal(h['indices'],a['ms']) and int(h['scale'])==S
assert F(str(h['error']))==sum(map(F,map(str,h['tile_errors'])),F(0))==F(meta['head_error'])
# Independent direct contraction: separate head and physical products,
# check actual stored floor and accumulate exact per-column squared error.
sq=[0]*n
def replay(Yh,Yp,Yout,accumulate):
    assert len(Yh)==len(Yp)==len(Yout)
    for k in range(0,len(Yout),128):
        A=Yh[k:k+128].astype(object)@q*S+Yp[k:k+128].astype(object)@x
        B=Yout[k:k+128].astype(object)
        R=A-B*S*S
        assert np.all(R>=0) and np.all(R<S*S)
        if accumulate:
            for j in range(n):sq[j]+=sum(int(v)**2 for v in R[:,j])
replay(pr['Y'][:,:119],pr['Y'][:,119:],old['yn'],False)
rows=[];common=F(0);total=F(0)
for chunk in ext['new_chunks']:
    lo,hi=chunk['lo'],chunk['hi'];v=np.load(eg/f'odd_U_{lo}_{hi}.npz');ids=v['indices']
    assert np.array_equal(ids,np.arange(lo,hi+1,2)) and int(v['scale'])==S
    assert v['Y'].shape==(len(ids),352)
    err=F(str(v['error']));assert err==sum(map(F,map(str,v['tile_errors'])),F(0))
    reserve=err-F(len(ids)*352,S);assert reserve>=0
    common+=reserve;total+=err;rows.extend(map(int,ids))
    ii=(ids-1)//2
    replay(h['Y'][ii],v['Y'],a['yn'][ii],True)
assert rows==list(range(16385,65536,2))
for j,v in enumerate(sq):
    diff=F(str(a['errors'][j]))-F(str(old['errors'][j]))
    assert diff>=0 and diff*diff>=F(v,S**6)
assert F(meta['physical_error'])==F(str(pr['U_error']))+total
# Independent attachment of block errors to the 65536 historical tail data.
oldbudget=read(c/'eb/generated/block-budget.json')
tail=np.load(dv/'merged/odd_merged_tail.npz');idx=list(map(int,tail['cutoff'])).index(65536)
assert F(meta['physical_tail'])==F(str(tail['tail'][idx]))
assert len(oldbudget['blocks'])==len(meta['physical_blocks'])==19
def su(v):
    assert v>=0
    g=10**30;z=isqrt(v.numerator*g*g//v.denominator)
    if F(z*z,g*g)<v:z+=1
    return F(z,g)
block=[]
for oldb,b in zip(oldbudget['blocks'],meta['physical_blocks']):
    assert (oldb['stage'],oldb['lo'],oldb['hi'])==(b['stage'],b['lo'],b['hi'])
    lo,hi=b['lo'],b['hi']
    assert F(b['finite_error'])==F(oldb['error'])+common+F(len(rows)*(hi-lo),S)
    tau=F(b['tail_squared']);assert tau>=0
    allow=F(b['finite_error'])+su(tau)
    assert allow==F(b['allowance']) and allow<=F(oldb['allowance'])/2
    block.append((lo,hi,allow,tau))
assert block[0][0]==0 and block[-1][1]==352
assert all(block[i][1]==block[i+1][0] for i in range(18))
assert sum(v[3] for v in block)==F(meta['physical_tail'])
aq=F(meta['head_error'])+su(F(meta['head_tail']))
ax=F(meta['physical_error'])+su(F(meta['physical_tail']))
assert aq==F(result['head_allowance']) and ax==F(result['physical_allowance'])
# Audit the retained circular row custody, not regenerate all primitive overlaps.
plan=read(eh/'circle-plan.json');ids=np.union1d(np.r_[a['js'],2*a['ms']],np.arange(1,32768,2))
assert len(ids)==plan['all_gram_rows'] and plan['parts']==(len(ids)+8191)//8192
parts=[np.load(eh/f'circle_part{k}.npz') for k in range(plan['parts'])]
assert np.array_equal(np.concatenate([v['ids'] for v in parts]),ids)
ix=np.flatnonzero((ids%2==1)&(ids<32768))
assert np.array_equal(circle['ids'],ids[ix])
for key in ('Z','Zerr'):assert np.array_equal(circle[key],np.concatenate([v[key] for v in parts])[ix]),key
assert digest(dz/'odd_joint_Uaction.npz')==plan['inverse_action_sha256']
er=su(sum(F(str(v))**2 for v in circle['reprerr']))
assert er==F(result['representation_error'])
# Independent exact assembly from frozen energy, new block allowances and
# published response matrices. H's underlying contractions remain inherited.
def cf(A,s):return [[F(int(v),s) for v in row] for row in A]
Q=cf(q.T@q,S*S);X=cf(x.T@x,S**4)
Xb=[cf(x[lo:hi].T@x[lo:hi],S**4) for lo,hi,_,_ in block]
# NPZ arbitrary-precision integer grids are serialized as decimal strings.
# astype(object) alone preserves strings; explicit Python int decoding is required.
obj=lambda A:np.vectorize(int,otypes=[object])(A)
en=np.load(dv/'merged/odd_BQ_energy.npz');K=obj(en['K']);D=obj(en['D']);KS=int(en['scale'])
assert K.shape==(352,119) and D.shape==(352,352) and KS>0
assert all(type(v) is int for v in K.flat) and all(type(v) is int for v in D.flat)
assert np.array_equal(D,D.T)
# Regression: exact integer conversion, multiplication, and contraction above int64.
test=obj(np.array([[str(10**80),'-3'],['7',str(-10**90)]]))
assert (test@np.array([[2],[5]],object)).tolist()==[[2*10**80-15],[14-5*10**90]]
print('EJ SYNTHESIS AND BUDGET REPLAY PASS; arbitrary precision energy grid decoded',flush=True)
E=cf((x.T@K@q+q.T@K.T@x)*S+x.T@D@x,KS*S**4)
eps=F(str(a['eps']));mu=F(str(a['mu']));ek=F(str(en['Kerror']));ed=F(str(en['Derror']))
assert mu>0 and ek>=0 and ed>=0
t0=a['t'].astype(object);alpha=a['a'].astype(object)
diag=cf(t0.T@((10**8-alpha)[:,None]*t0),10**8*S*S)
for i in range(n):
    for j in range(n):E[i][j]+=(1-eps)*Q[i][j]
qd=max(sum(abs(Q[i][j]-(i==j)) for j in range(n)) for i in range(n))
def mat(r,k):
    v=[[F(x) for x in row] for row in r[k]]
    assert len(v)==n and all(len(row)==n for row in v)
    assert all(v[i][j]==v[j][i] for i in range(n) for j in range(n))
    return v
def bareiss(A):
    A=[row[:] for row in A];prev=1;ds=[]
    for k in range(n):
        d=A[k][k];assert d>0;ds.append(d)
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j],rem=divmod(d*A[i][j]-A[i][k]*A[k][j],prev);assert rem==0
        prev=d
    return ds
omit={'matrix','H','R','P','F_lower','Gu_upper','Bu_upper','Bu_lower','Cross_upper'}
checks=[]
paths=sorted([*ei.glob('matched-eb-parameters-center-*.json'),*ei.glob('optimized-*.json')])
assert len(paths)==len(result['comparisons'])
published=[]
for path in paths:
    r=read(path);assert r['n']==n and r['mode']==1 and r['split']
    p={k:v for k,v in r.items() if k not in omit};published.append(p)
    assert p in result['comparisons']
    th,de,t,eta=map(F,r['parameters']);w=list(map(F,r['weights']));center=F(r['center'])
    assert min(th,de,t,eta,*w)>0 and sum(w)==1 and len(w)==21 and 0<=center<=F(1,2)
    assert qd==F(r['qdefect'])<1
    M,H,R,P,V,G,B,C=(mat(r,k) for k in ('matrix','H','R','P','F_lower','Gu_upper','Bu_lower','Cross_upper'))
    for i in range(n):
        for j in range(n):
            assert R[i][j]==H[i][j]-diag[i][j]-V[i][j]
            Ac=(1-2*center)*G[i][j]+center*C[i][j]-2*center*B[i][j]+center*center*R[i][j]
            assert P[i][j]==(1+t)*Ac+(1+1/t)*(1-center)**2*R[i][j]
            DD=aq*aq/w[0]*Q[i][j]+sum(b[2]**2/w[k+1]*Xb[k][i][j] for k,b in enumerate(block))+er*er/w[-1]*(i==j)
            EE=E[i][j]-ek*(eta*Q[i][j]+X[i][j]/eta)-ed*X[i][j]
            assert M[i][j]==EE-(H[i][j]+de*P[i][j]+(1+1/de)*DD)/mu
    if not r['positive']:continue
    lower=F(r['lower']);assert lower>0
    if 'margin' not in r:continue
    margin=F(r['margin']);assert 0<margin<lower
    shift=margin*(1+qd)
    T=[[(M[i][j]-shift*(i==j))*S//1-n*(i==j) for j in range(n)] for i in range(n)]
    assert [T[i][j] for i in range(n) for j in range(i,n)]==r['upper_triangle']
    ds=bareiss(T)
    checks.append(dict(file=path.name,lower=str(lower),margin=str(margin),leading_minors=list(map(str,ds))))
assert max(F(v['lower']) for v in published)==F(result['best']['lower'])
assert result['best'] in published
accepted=bool(result['best']['positive'] and any(F(v['lower'])==F(result['best']['lower']) for v in checks))
report=dict(scope='independent synthesis contraction, directed-budget attachment and terminal rational matrix/integer positivity audit; primitive circle and analytical enclosures not independently regenerated',full_synthesis_rows=32768,columns=n,physical_blocks=19,comparison_assemblies=len(paths),positive_checks=checks,leading_minor_count=n*len(checks),best_margin_audited=accepted)
assert accepted and len(checks)==8 and n*len(checks)==248
(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
print('EJ AUDIT',len(paths),len(checks),accepted,flush=True)
