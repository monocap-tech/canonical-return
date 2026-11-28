"""Combine recovered and new shards; no common-proxy assumption in the primary bound."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,sys
root,logs,prior=map(Path,sys.argv[1:4])
old=json.loads((root/'cr_alp_058'/'verification_result.json').read_text())['2']
second=json.loads((prior/'cr_alp_059/evidence/second_result.json').read_text())['2']
third=json.loads((prior/'cr_alp_061/evidence/third_result.json').read_text())['2']
fourth=json.loads((prior/'cr_alp_063/evidence/fourth_result.json').read_text())['2']
target=json.loads((prior/'cr_alp_064/fifth_moment_targets.json').read_text())
assert target['sector']=='2'
ell,u,g,s=map(F,('0.255','5.505','2.88','2.625'))
def parse(path):
    rows={}
    for line in path.read_text().splitlines():
        x=line.split()
        if len(x)==5 and x[-2:]==['/','1e15']:
            assert x[0] not in rows
            rows[x[0]]=tuple(F(int(v),10**15) for v in x[1:3])
    assert all(l<=h for l,h in rows.values())
    return rows
base=prior/'cr_alp_065/evidence/fifth_2_partial_timeout.log'
rows=parse(base)
prefix=set(rows)-{'fifth_corr0','fifth_corr1'}
assert len(prefix)==11 and set(rows)==prefix|{'fifth_corr0','fifth_corr1'}
fingerprints=[]
for i in range(2,7):
    path=logs/f'fifth_shift_{i}.log'
    shard=parse(path)
    key=f'fifth_corr{i}'
    assert set(shard)==prefix|{key}
    assert all(shard[k]==rows[k] for k in prefix)
    exact=[line for line in path.read_text().splitlines() if line.startswith('midpoints65 ')]
    assert len(exact)==1 and len(exact[0].split())>3
    fingerprints.append(exact[0])
    rows[key]=shard[key]
assert len(set(fingerprints))==1
# This fingerprint matches the new projection coordinates, not the old timeout job.
# Even identical rounded norm rows cannot prove equality of all midpoint functions.

priorrows=parse(prior/'cr_alp_063/evidence/fourth_2.log')
assert len(priorrows)==10
assert all(rows[k]==v for k,v in priorrows.items())
assert set(rows)==set(priorrows)|{'projection_error2',*(f'fifth_corr{i}' for i in range(7))}
def overlap(x,y):return max(x[0],y[0])<=min(x[1],y[1])
lines=(root/'cr_alp_058/mesh.txt').read_text().splitlines()
n=int(lines[0]);weights=[tuple(map(F,l.split()[2:])) for l in lines[n+1:n+8]]
assert len(weights)==7 and all(0<l<=h for l,h in weights)
assert 2*sum(h for l,h in weights)<8
# A v0 minus finite correction: amplitude <1e20, derivative <1e23,
# fewer than 16384 jump neighborhoods. Paired shift perturbation.
assert 8*(2*10**23*F(1,10**70)+2*10**20*F(128,10**35))<F(1,10**11)
eta2=rows['projection_error2'][1]
assert rows['projection_error2'][0]<=0<=eta2<F(1,10)
B=F(fourth['proxy_norm_bound'])+eta2
assert B<F(11,10)
delta=F(fourth['second_application_error'])+eta2

corr=[[-2*w*v for w,v in product(weights[i],rows[f'fifth_corr{i}'])] for i in range(7)]
recon=(sum(map(min,corr)),sum(map(max,corr)))
p=recon
# Each correlation may use its own w0, independently within delta of p^2 r.
# Each one-sided translation is a contraction. The weighted displacement
# error bound applies termwise, with B<1.1; no shared proxy needed.
weight_sum=2*sum(h for l,h in weights)
robust_err=weight_sum*(2*B*delta+delta**2)+F(2,10**11)
# Alternative ONLY if all seven runs have the identical exact physical proxy.
common_err=s*(2*B*delta+delta**2)+F(2,10**11)
support=tuple(map(F,target['support_outer']))
def enclosure(err):
    raw=(p[0]-err,p[1]+err)
    assert overlap(raw,support)
    return (max(raw[0],support[0]),min(raw[1],support[1]))

G,M,K=(tuple(map(F,old[x])) for x in ('true_G','true_M','K'))
N=tuple(map(F,second['true_N']));T=tuple(map(F,third['true_T']));S=tuple(map(F,fourth['true_S']))
x,y=map(F,target['trial_upper']);a,b,c=map(F,target['trial_lower'])
def upper(G,M,N,T,S,V):
    h=s*s*G-N;j=s*s*M-T;k=s*s*N-S
    return (g*G-M-2*(h*x+j*y)+(g*h+j)*x*x+2*(g*j+k)*x*y+(g*k+s*s*T-V)*y*y)/(ell*u)
def lower(G,M,N,T,S,V):
    m=[G,M,N,T,S,V];z=[a,b,c]
    return 2*(G*a+M*b+N*c)-sum(z[i]*z[j]*(g*m[i+j]+m[i+j+1]) for i in range(3) for j in range(3))
def decide(err):
    F5=enclosure(err)
    corners=list(product(G,M,N,T,S,F5))
    hi=min(max(upper(*v) for v in corners),F(fourth['inverse_bounds'][1]))
    lo=max(min(lower(*v) for v in corners),F(fourth['inverse_bounds'][0]))
    assert lo<=hi
    status='PASS_WITNESS' if hi<K[0] else 'FAIL_COMPARISON_WITNESS' if lo>K[1] else 'UNDECIDED'
    return {'transport':str(err),'true_F':list(map(str,F5)),'inverse_bounds':list(map(str,(lo,hi))),'status':status}
out={'sector':'2','weighted_proxy_correlation_sum':list(map(str,p)),'second_projection_error':str(eta2),'proxy_norm_bound':str(B),'total_L2_error':str(delta),'paired_weight_sum_upper':str(weight_sum),'K':list(map(str,K)),'new_shard_midpoint_coordinates_match':True,'primary_independent_proxy_transport':decide(robust_err),'conditional_identical_proxy_transport':decide(common_err)}
for key in ('primary_independent_proxy_transport','conditional_identical_proxy_transport'):
    v=out[key]
    print(key,'F',*[float(F(x)) for x in v['true_F']],'inverse',*[float(F(x)) for x in v['inverse_bounds']],v['status'])
(logs/'fifth_sharded_result.json').write_text(json.dumps(out,indent=2)+'\n')
