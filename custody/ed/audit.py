from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
root=Path(__file__).resolve().parents[2];src=root/'custody/ec/generated';out=root/'custody/ed/generated';out.mkdir(parents=True,exist_ok=True)
d=json.loads((src/'family-obstruction.json').read_text())
v=json.loads((src/'centered-control.json').read_text())
z=[F(x,d['witness_scale']) for x in d['witness']]
n=len(z);assert n==31 and d['witness_scale']==10**12
N=sum(x*x for x in z);assert N==F(d['witness_norm_squared'])>0
c={k:F(x) if isinstance(x,str) else list(map(F,x)) for k,x in d['ray_constants'].items()}
def quad(key):
    A=[[F(x) for x in row] for row in v[key]]
    assert len(A)==n and all(len(row)==n for row in A)
    assert all(A[i][j]==A[j][i] for i in range(n) for j in range(n))
    return sum(z[i]*A[i][j]*z[j] for i in range(n) for j in range(n))
th,de,t,eta=map(F,v['parameters']);w=list(map(F,v['weights']))
assert min(th,de,t,eta,*w)>0 and sum(w)==1 and len(w)==21
assert v['center']=='1/2' and v['lower']=='-3/5000'
h=c['base']+c['r']*(th*c['W']+c['U']/th)
rr=h-c['diag']-c['F_lower']
pp=(1+t)*(c['Cross']/2-c['Bu_lower']+rr/4)+(1+1/t)*rr/4
dd=sum(x/y for x,y in zip(c['block_terms'],w))
ee=c['energy']-c['eK']*(eta*c['Q']+c['X']/eta)-c['eD']*c['X']
assert h==quad('H') and rr==quad('R') and pp==quad('P')
for k,key in [('F_lower','F_lower'),('Gu','Gu_upper'),('Cross','Cross_upper'),('Bu_lower','Bu_lower')]:
    assert c[k]==quad(key),k
L=quad('matrix');assert L==ee-(h+de*pp+(1+1/de)*dd)/c['mu'] and L<0
# Independent lower square-root implementation: floor on a finer dyadic grid.
S=2**100
def lo(x):
    assert x>=0
    k=isqrt(x.numerator*S*S//x.denominator)
    assert k*k*x.denominator<=x.numerator*S*S<(k+1)**2*x.denominator
    return F(k,S)
assert min(c[k] for k in ['Q','X','W','U','r','eK','eD','Gu'])>=0 and c['mu']>0
assert min(c['scalar_terms']+c['block_terms'])>=0
assert len(c['scalar_terms'])==3 and len(c['block_terms'])==21
H=c['base']+2*c['r']*lo(c['W']*c['U'])
R=max(F(0),H-c['diag']-c['F_lower'])
E=c['energy']-2*c['eK']*lo(c['Q']*c['X'])-c['eD']*c['X']
# Different cover: linear and quadratic terms separately bounded on each cell.
# This avoids reuse of the producer's vertex calculation.
slope=c['Cross']-2*c['Bu_lower']-2*c['Gu']
P1=None
for i in range(2048):
    a=F(i,4096);b=F(i+1,4096)
    A=max(F(0),c['Gu']+min(slope*a,slope*b)+R*a*a)
    p=(lo(A)+(1-b)*lo(R))**2
    P1=p if P1 is None else min(P1,p)
P0=max(F(0),H-c['diag'])
def bound(terms,P):
    D=sum(lo(x) for x in terms)**2
    return (E-(H+D+2*lo(P*D))/c['mu'])/N
checks=[]
for family in ['scalar','block']:
    for mode,P in [('mode0',P0),('centered-mode1',P1)]:
        x=bound(c[family+'_terms'],P)
        assert x<0
        checks.append(dict(family=family,mode=mode,upper=str(x),decimal=float(x)))
# Sensitivity of the necessary-condition witness, not an achievable margin.
rows=[]
terms=c['block_terms']
for head in [F(1),F(3,4),F(1,2),F(1,4),F(0)]:
    for phys in [F(1),F(3,4),F(1,2),F(1,4),F(0)]:
        new=[head*head*terms[0]]+[phys*phys*x for x in terms[1:-1]]+[terms[-1]]
        x=bound(new,P1)
        rows.append(dict(head_norm_factor=str(head),physical_norm_factor=str(phys),upper=str(x),decimal=float(x),still_excluded=bool(x<0)))
report=dict(scope='independent terminal witness arithmetic on stored ray constants; analytical enclosures and all raw ray constants not independently regenerated',witness_matrix_reconstruction=True,interval_cells=2048,checks=checks,all_excluded=True,sensitivity=rows)
(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
print('ED AUDIT',[(x['family'],x['mode'],x['decimal']) for x in checks],flush=True)
print('ED SENSITIVITY',[(x['head_norm_factor'],x['physical_norm_factor'],x['decimal']) for x in rows],flush=True)
