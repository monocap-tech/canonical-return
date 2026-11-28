"""CR-ALP-063 exact rational fourth-moment transport and frozen trial checks."""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import json,sys
root,logs,prior=map(Path,sys.argv[1:4])
old=json.loads((root/'cr_alp_058'/'verification_result.json').read_text())
second=json.loads((prior/'cr_alp_059'/'evidence'/'second_result.json').read_text())
third=json.loads((prior/'cr_alp_061'/'evidence'/'third_result.json').read_text())
targets=json.loads((prior/'cr_alp_061'/'evidence'/'fourth_moment_trials.json').read_text())
support=json.loads((prior/'cr_alp_061'/'evidence'/'fourth_moment_feasibility.json').read_text())
ell,u,g,s=map(F,('0.255','5.505','2.88','2.625'))
assert ell==g-s and u==g+s
# Independently verify displacement allowance for the frozen midpoint projection v0.
assert 8*(2*10**22*F(1,10**70)+2*10**19*F(32,10**35))<F(1,10**12)
lines=(root/'cr_alp_058'/'mesh.txt').read_text().splitlines()
n=int(lines[0]);weights=[tuple(map(F,l.split()[2:])) for l in lines[n+1:n+8]]
assert len(weights)==7 and all(0<l<=h for l,h in weights)
assert 2*sum(h for l,h in weights)<8
def overlap(x,y):return max(x[0],y[0])<=min(x[1],y[1])
def sqrt_upper(x):
    assert x>=0
    D=10**30
    y=x*D**2
    return F(isqrt(y.numerator//y.denominator)+1,D)
result={}
for k in ('1','2'):
    rows={}
    for line in (logs/f'fourth_{k}.log').read_text().splitlines():
        x=line.split()
        if len(x)==5 and x[-2:]==['/','1e15']:
            assert x[0] not in rows
            rows[x[0]]=tuple(F(int(v),10**15) for v in x[1:3])
    assert set(rows)=={'proxy_G','V2','sine2','bridge2','N_proxy','projection_error','W2','sine4','bridge4','S_proxy'}
    assert all(l<=h for l,h in rows.values())
    assert 0<=rows['N_proxy'][0]<=rows['N_proxy'][1]<F(11,20)**2
    assert overlap(rows['N_proxy'],tuple(map(F,second[k]['proxy_N'])))
    eta=rows['projection_error'][1]
    assert rows['projection_error'][0]<=0<=eta<F(1,10)
    # The same frozen coordinate selection is regenerated in the same order.
    assert eta==F(third[k]['projection_error'])
    w,a,b,p=map(rows.get,('W2','sine4','bridge4','S_proxy'))
    assert 0<=p[0]<=p[1]<10 and overlap(p,(max(F(0),w[0]-a[1]-b[1]),w[1]-a[0]-b[0]))
    eps=F(old[k]['proxy_error_upper'])
    delta=s*eps+F(1,10**13)+eta
    delta2=s*delta+F(1,10**12)
    # ||Q P Q P r-Q A v0||<=delta2. Use the observed norm upper bound.
    bound=sqrt_upper(p[1])
    transport=2*bound*delta2+delta2**2
    S=(max(F(0),p[0]-transport),p[1]+transport)
    admissible=(F(support[k]['support_lower_boundary'][0]),F(support[k]['support_upper_boundary'][1]))
    assert overlap(S,admissible)
    S=(max(S[0],admissible[0]),min(S[1],admissible[1]))
    G,M,K=(tuple(map(F,old[k][x])) for x in ('true_G','true_M','K'))
    N=tuple(map(F,second[k]['true_N']))
    T=tuple(map(F,third[k]['true_T']))
    y,beta=map(F,targets[k]['trial_upper']); z,delta_trial=map(F,targets[k]['trial_lower'])
    def upper(gg,mm,nn,tt,ss):
        b0=s*gg+mm;b1=s*mm+nn
        A=g*s*gg+u*mm+nn;B=g*s*mm+u*nn+tt;C=g*s*nn+u*tt+ss
        return (gg-2*(b0*y+b1*beta)+A*y*y+2*B*y*beta+C*beta*beta)/ell
    def lower(gg,mm,nn,tt,ss):
        d0=s*gg-mm;d1=s*mm-nn
        D=g*s*gg-ell*mm-nn;E=g*s*mm-ell*nn-tt;C=g*s*nn-ell*tt-ss
        return (gg+2*(d0*z+d1*delta_trial)-D*z*z-2*E*z*delta_trial-C*delta_trial*delta_trial)/u
    corners=list(product(G,M,N,T,S))
    hi=min(max(upper(*x) for x in corners),F(third[k]['inverse_bounds'][1]))
    lo=max(min(lower(*x) for x in corners),F(third[k]['inverse_bounds'][0]))
    assert lo<=hi
    status='PASS_WITNESS' if hi<K[0] else 'FAIL_COMPARISON_WITNESS' if lo>K[1] else 'UNDECIDED'
    result[k]={'S_proxy':list(map(str,p)),'projection_error':str(eta),'second_application_error':str(delta2),'proxy_norm_bound':str(bound),'transport':str(transport),'true_S':list(map(str,S)),'inverse_bounds':list(map(str,(lo,hi))),'K':list(map(str,K)),'status':status}
    print(k,'S',*[float(x) for x in S],'inverse',float(lo),float(hi),'K',*[float(x) for x in K],status)
(logs/'fourth_result.json').write_text(json.dumps(result,indent=2)+'\n')
