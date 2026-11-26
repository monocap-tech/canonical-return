from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,json
root=Path(__file__).resolve().parents[2];out=root/'custody/eg/generated';old=np.load(root/'custody/dv/generated/merged/odd_combined_projection.npz')
bp=np.load(root/'ce-output/pristine/cr_int_t1d_bp/odd_head_full.npz');h=np.load(out/'odd_head_full.npz')
assert np.array_equal(h['indices'],np.arange(1,65536,2))
assert h['Y'].shape==(32768,119) and int(h['scale'])==int(bp['scale'])==10**12
assert np.array_equal(h['Y'][:8192],bp['Y']) and np.array_equal(bp['Y'],old['Y'][:,:119])
eh=F(str(h['error']));assert eh==sum(map(lambda x:F(str(x)),h['tile_errors']),F(0))
assert F(str(h['tile_errors'][0]))==F(str(bp['error']))==F(str(old['head_error']))
pairs=[(16385,18431),(18433,26623),(26625,34815),(34817,43007),(43009,51199),(51201,59391),(59393,65023),(65025,65535)]
rows=[];errors=[];blocks=json.loads((root/'custody/ee/generated/census.json').read_text())
budget=json.loads((root/'custody/eb/generated/block-budget.json').read_text())
assert len(budget['blocks'])==19
S=10**12
# The inherited tile bound includes a deterministic entrywise rounding allowance
# of rows*352/S; restrict this term to block columns, retaining the entire
# remaining non-rounding allowance for each block.
common=F(0);totalrows=0
for lo,hi in pairs:
    a=np.load(out/f'odd_U_{lo}_{hi}.npz');ids=a['indices'];Y=a['Y'];err=F(str(a['error']))
    assert np.array_equal(ids,np.arange(lo,hi+1,2)) and Y.shape==(len(ids),352) and int(a['scale'])==S
    assert err==sum(map(lambda x:F(str(x)),a['tile_errors']),F(0))
    extra=err-F(len(ids)*352,S);assert extra>=0
    common+=extra;totalrows+=len(ids);rows.extend(ids.tolist());errors.append(err)
assert rows==list(range(16385,65536,2)) and totalrows==24576
physical=F(str(old['U_error']))+sum(errors)
T=np.load(root/'custody/dv/generated/merged/odd_merged_tail.npz')
j=list(map(int,T['cutoff'])).index(65536)
headtail=np.load(root/'ce-output/pristine/cr_int_t1d_bo/odd_tail.npz')
k=list(map(int,headtail['cutoff'])).index(65536)
def su(x):
    s=10**30;n=-((-x.numerator*s*s)//x.denominator);r=isqrt(n);a=F(r+(r*r<n),s);assert a*a>=x;return a
head_allow=eh+su(F(str(headtail['tail'][k])))
head_old=F(blocks['head']['old_allowance'])
summary=[]
for b in budget['blocks']:
    width=b['hi']-b['lo']
    rec=next(x for x in blocks['blocks'] if x['stage']==b['stage'] and x['cutoff']==65536)
    e=F(b['error'])+common+F(totalrows*width,S)
    tail=F(rec['tail_squared']);allow=e+su(tail);old_a=F(b['allowance'])
    summary.append(dict(stage=b['stage'],lo=b['lo'],hi=b['hi'],finite_error=str(e),tail_squared=str(tail),allowance=str(allow),ratio=float(allow/old_a),half_target_met=allow<=old_a/2))
report=dict(scope='directed physical and normalized-head finite action and inherited tails; physical synthesis and circle symbol not yet regenerated',cutoff=65536,head_rows=len(h['Y']),physical_new_rows=totalrows,physical_columns=352,head_columns=119,head_error=str(eh),head_tail=str(F(str(headtail['tail'][k]))),head_allowance=str(head_allow),head_ratio=float(head_allow/head_old),head_half_target_met=head_allow<=head_old/2,physical_finite_error=str(physical),physical_tail=str(F(str(T['tail'][j]))),physical_allowance=str(physical+su(F(str(T['tail'][j])))),physical_blocks=summary,all_block_half_targets_met=all(b['half_target_met'] for b in summary),new_chunks=[dict(lo=l,hi=h) for l,h in pairs])
(out/'extension.json').write_text(json.dumps(report,indent=2)+'\n')
print('EG COMPLETE COVERAGE',totalrows,'head ratio',report['head_ratio'],'maximum block ratio',max(v['ratio'] for v in summary),flush=True)
