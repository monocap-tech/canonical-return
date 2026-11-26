"""EN: attach one new 16-column block to the frozen extended action.
No positivity assertion; new directed action and inherited analytic tail only.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,json,hashlib
root=Path(__file__).resolve().parents[2];c=root/'custody';out=c/'en/generated';new=out/'new';merged=out/'merged';eg=c/'eg/generated'
def auth(p):
    for line in (p/'sha256.txt').read_text().splitlines():
        digest,name=line.split(maxsplit=1)
        with (p/name.removeprefix('./')).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==digest
for stage in ('dv','eg','em'):auth(c/stage/'generated')
em=c/'em/generated';receiver=json.loads((em/'receiver-audit.json').read_text())
assert receiver['columns']==16 and F(receiver['rho'])<F(1,100)
assert np.array_equal(np.load(em/'odd_new.npy'),receiver['new_indices'])
old=json.loads((eg/'extension.json').read_text())
assert old['cutoff']==65536 and old['physical_columns']==352 and old['all_block_half_targets_met']
S=10**12
def su(v):
    assert v>=0
    g=10**30;z=isqrt(v.numerator*g*g//v.denominator)
    if F(z*z,g*g)<v:z+=1
    return F(z,g)
pairs=[(1,2047),(2049,6143),(6145,11263),(11265,16383)]+[(lo,min(lo+8190,65535)) for lo in range(16385,65536,8192)]
ids=[];ys=[];err=F(0);prefixerr=F(0)
for lo,hi in pairs:
    a=np.load(new/f'odd_U_{lo}_{hi}.npz')
    assert np.array_equal(a['indices'],np.arange(lo,hi+1,2)) and a['Y'].shape==(len(a['indices']),16) and int(a['scale'])==S
    e=F(str(a['error']));assert e==sum(map(F,map(str,a['tile_errors'])),F(0)) and e>=F(a['Y'].size,S)
    err+=e
    if hi<16384:prefixerr+=e
    ids.extend(map(int,a['indices']));ys.append(a['Y'])
assert ids==list(range(1,65536,2))
Y=np.concatenate(ys)
prior=np.load(c/'dv/generated/merged/odd_combined_projection.npz');p=np.load(merged/'odd_combined_projection.npz')
assert p['Y'].shape==(8192,487)
assert np.array_equal(p['Y'][:,:471],prior['Y']) and np.array_equal(p['Y'][:,471:],Y[:8192])
assert F(str(p['U_error']))==F(str(prior['U_error']))+prefixerr
h=np.load(eg/'odd_head_full.npz');assert np.array_equal(h['Y'][:8192],prior['Y'][:,:119])
# Preserve all old extended chunks byte-for-byte; new block is stored separately.
blocks=old['physical_blocks'];assert len(blocks)==19 and blocks[-1]['hi']==352
tail=np.load(em/'odd_U_tail.npz');i=list(map(int,tail['cutoff'])).index(65536)
tau=F(str(tail['tail'][i]));t=np.load(merged/'odd_merged_tail.npz');j=list(map(int,t['cutoff'])).index(65536)
assert F(str(t['tail'][j]))==F(old['physical_tail'])+tau
allow=err+su(tau)
blocks=blocks+[dict(stage='en',lo=352,hi=368,finite_error=str(err),tail_squared=str(tau),allowance=str(allow))]
assert sum(F(b['tail_squared']) for b in blocks)==F(str(t['tail'][j]))
co=np.load(merged/'odd_merged_corrector.npz');en=np.load(merged/'odd_BQ_energy.npz')
assert co['coefficient'].shape==(507,368) and en['K'].shape==(368,119) and en['D'].shape==(368,368)
# Preserve the old energy submatrices within the sum of both interval errors.
a=np.load(c/'dv/generated/merged/odd_BQ_energy.npz')
def norm(A):
    v=sum(int(x)**2 for x in A.flat);z=isqrt(v);return z+(z*z<v)
obj=lambda A:np.vectorize(int,otypes=[object])(A)
s0,s1=int(a['scale']),int(en['scale'])
for key,sl,ek in [('K',(slice(0,352),slice(None)),'Kerror'),('D',(slice(0,352),slice(0,352)),'Derror')]:
    diff=obj(en[key])[sl]*s0-obj(a[key])*s1
    assert F(norm(diff),s0*s1)<=F(str(a[ek]))+F(str(en[ek])),key
report=dict(scope='new directed physical block and merged energy; not terminal positivity',cutoff=65536,physical_columns=368,physical_support=507,new_columns=16,new_rows=len(ids),new_finite_error=str(err),new_tail_squared=str(tau),new_allowance=str(allow),aggregate_finite_error=str(F(old['physical_finite_error'])+err),aggregate_tail_squared=str(F(old['physical_tail'])+tau),head_error=old['head_error'],head_tail=old['head_tail'],physical_blocks=blocks,old_energy_overlap=True,old_action_reused=True,new_chunks=[dict(lo=lo,hi=hi) for lo,hi in pairs])
(out/'extension.json').write_text(json.dumps(report,indent=2)+'\n')
print('EN EXTENDED BLOCK PASS',len(ids),len(blocks),float(allow),flush=True)
