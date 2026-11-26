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
# No new directed action is performed: census only.
from math import isqrt
import shutil
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=workspace.parent.parent;out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
d=loadjoint()
old=json.loads((root/'custody/eb/generated/block-budget.json').read_text())
assert old['cutoff']=='16384' and len(old['blocks'])==19
def f(x):return F(str(x))
def load(s,p):return np.load(root/'custody'/s/'generated'/p)
S0=10**24
def su(x):
    assert x>=0
    k=isqrt(x.numerator*S0*S0//x.denominator)
    a=F(k,S0)
    if a*a<x:a+=F(1,S0)
    assert a*a>=x
    return a
stages=('cm','co','cq','cr','ct','cv','cy','da','db','dc','de','dg','di','dk','dn','ds','du','dv')
full=load('dv','merged/odd_merged_tail.npz');cut=list(map(int,full['cutoff']))
assert cut[0]==16384 and 65536 in cut
T=[list(map(f,load(s,'new/odd_U_tail.npz')['tail'])) for s in stages]
for s in stages:assert list(map(int,load(s,'new/odd_U_tail.npz')['cutoff']))==cut
first=list(map(f,load('cm','merged/odd_merged_tail.npz')['tail']))
base=[x-y for x,y in zip(first,T[0])]
alltail=[base]+T
assert all(sum(v[k] for v in alltail)==f(full['tail'][k]) for k in range(len(cut)))
assert all(v>=0 for row in alltail for v in row)
rows=[]
for ob,tau in zip(old['blocks'],alltail):
    e=f(ob['error']);aold=f(ob['allowance'])
    assert tau[0]==f(ob['tail']) and e>=0 and aold>=e and (aold-e)**2>=tau[0]
    for k,N in enumerate(cut):
        tailnorm=su(tau[k]);target=aold/2
        rows.append(dict(stage=ob['stage'],cutoff=N,old_allowance=str(aold),old_finite_error=str(e),tail_squared=str(tau[k]),tail_norm_upper=str(tailnorm),tail_only_ratio=float(tailnorm/aold),hypothetical_same_finite_error_ratio=float((e+tailnorm)/aold),new_finite_error_cap_for_half=str(target-tailnorm),cap_positive=target>=tailnorm))
# Identify actual head error/tail data rather than assuming its decay matches the correctors.
head=dict(old_eY=str(d['eY']),old_tau=str(d['tau']),old_allowance=str(d['eY']+sqrtup(d['tau'])),target_half=str((d['eY']+sqrtup(d['tau']))/2))
assert d['eY']>=0 and d['tau']>=0
assert f(old['scalar_allowance'])==d['eU']+sqrtup(d['tauU'])
# Preserve source files and exact lines relevant to changing cutoffs.
srcdir=out/'sources';srcdir.mkdir(exist_ok=True)
inventory=[]
wanted={'bq_action.py','bq_tail.py','combine.py','bq_project.py','bp_action.py','bp_tail.py','bp_project.py','bs_trial.py','bw_prepare.py','bs_circle.py','bs_combine.py','bv_joint_action.py','bv_joint_trial.py','bv_response.py','cb_matrix.py'}
for p in sorted(workspace.rglob('*.py')):
    if p.name not in wanted:continue
    rel=p.relative_to(workspace);dest=srcdir/rel;dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(p,dest)
    lines=p.read_text().splitlines()
    inventory.append(dict(path=str(rel),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),cutoff_lines=[dict(line=i+1,text=l) for i,l in enumerate(lines) if any(k in l for k in ('16384','8192','32768','65536','tauU','eY','head_error','tau','cutoff'))]))
headfiles=[]
for p in sorted(workspace.rglob('*.npz')):
    if 'tail' not in p.name.lower():continue
    try:
        with np.load(p) as v:
            entry=dict(path=str(p.relative_to(workspace)),keys=list(v.files))
            for k in ('cutoff','tail','tau','error'):
                if k in v.files and v[k].size<=32:entry[k]=list(map(str,v[k].reshape(-1)))
            headfiles.append(entry)
    except (ValueError, OSError) as e:
        headfiles.append(dict(path=str(p.relative_to(workspace)),error=str(e)))
report=dict(scope='authenticated historical tail bounds and hypothetical finite error budgets only; no regenerated higher cutoff action or certificate',head=head,cutoffs=cut,blocks=rows,source_inventory=inventory,tail_file_inventory=headfiles)
(out/'census.json').write_text(json.dumps(report,indent=2)+'\n')
print('EE HEAD',head,flush=True)
print('EE PHYSICAL 65536',[(v['stage'],v['tail_only_ratio'],v['hypothetical_same_finite_error_ratio'],v['cap_positive']) for v in rows if v['cutoff']==65536],flush=True)
print('EE SOURCE INVENTORY',len(inventory),'TAIL FILES',len(headfiles),flush=True)
