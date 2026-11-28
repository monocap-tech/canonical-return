"""Exact even-only fifth-moment transport and frozen weighted inverse tests."""
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
rows=parse(logs/'fifth_2.log')
priorrows=parse(prior/'cr_alp_063/evidence/fourth_2.log')
assert len(priorrows)==10
assert all(rows[k]==v for k,v in priorrows.items())
assert set(rows)==set(priorrows)|{'projection_error2','F_proxy',*(f'fifth_corr{i}' for i in range(7))}
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
err=s*(2*B*delta+delta**2)+F(2,10**11)
p=rows['F_proxy']
corr=[[-2*w*v for w,v in product(weights[i],rows[f'fifth_corr{i}'])] for i in range(7)]
recon=(sum(map(min,corr)),sum(map(max,corr)))
assert overlap(p,recon)
raw=(p[0]-err,p[1]+err)
support=tuple(map(F,target['support_outer']))
assert overlap(raw,support)
F5=(max(raw[0],support[0]),min(raw[1],support[1]))
G,M,K=(tuple(map(F,old[x])) for x in ('true_G','true_M','K'))
N=tuple(map(F,second['true_N']));T=tuple(map(F,third['true_T']));S=tuple(map(F,fourth['true_S']))
x,y=map(F,target['trial_upper']);a,b,c=map(F,target['trial_lower'])
def upper(G,M,N,T,S,V):
    h=s*s*G-N;j=s*s*M-T;k=s*s*N-S
    return (g*G-M-2*(h*x+j*y)+(g*h+j)*x*x+2*(g*j+k)*x*y+(g*k+s*s*T-V)*y*y)/(ell*u)
def lower(G,M,N,T,S,V):
    m=[G,M,N,T,S,V];z=[a,b,c]
    return 2*(G*a+M*b+N*c)-sum(z[i]*z[j]*(g*m[i+j]+m[i+j+1]) for i in range(3) for j in range(3))
corners=list(product(G,M,N,T,S,F5))
hi=min(max(upper(*v) for v in corners),F(fourth['inverse_bounds'][1]))
lo=max(min(lower(*v) for v in corners),F(fourth['inverse_bounds'][0]))
assert lo<=hi
status='PASS_WITNESS' if hi<K[0] else 'FAIL_COMPARISON_WITNESS' if lo>K[1] else 'UNDECIDED'
out={'sector':'2','F_proxy':list(map(str,p)),'second_projection_error':str(eta2),'proxy_norm_bound':str(B),'total_L2_error':str(delta),'transport':str(err),'true_F':list(map(str,F5)),'inverse_bounds':list(map(str,(lo,hi))),'K':list(map(str,K)),'status':status}
print('2 F',*[float(v) for v in F5],'inverse',float(lo),float(hi),'K',*[float(v) for v in K],status)
(logs/'fifth_result.json').write_text(json.dumps(out,indent=2)+'\n')
