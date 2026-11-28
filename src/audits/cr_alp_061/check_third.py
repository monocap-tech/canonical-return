"""Exact transport and cubic/variational witness tests; no Lean claim."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import json,sys,runpy
root=Path(sys.argv[1]); logs=Path(sys.argv[2]); old=json.loads((root/'cr_alp_058'/'verification_result.json').read_text())
second=json.loads(Path(sys.argv[3]).read_text())
targets=json.loads(Path(sys.argv[4]).read_text())
ell,u,g,s=map(F,('0.255','5.505','2.88','2.625'))
# Parse exact rational weight enclosures from frozen mesh.
lines=(root/'cr_alp_058'/'mesh.txt').read_text().splitlines()
n=int(lines[0]);weights=[tuple(map(F,l.split()[2:])) for l in lines[n+1:n+8]]
assert len(weights)==7 and all(0<l<=h for l,h in weights)
assert 2*sum(h for l,h in weights)<8
# v0 amplitude <1e19 and derivative <1e22; <=1024 jump neighborhoods.
# Shift error in L2 with the sum of paired weights <8.
assert 8*(2*10**22*F(1,10**70)+2*10**19*F(32,10**35))<F(1,10**12)
result={}
for k in ('1','2'):
    rows={}
    for line in (logs/f'third_{k}.log').read_text().splitlines():
        x=line.split()
        if len(x)==5 and x[-2:]==['/','1e15']:
            assert x[0] not in rows
            rows[x[0]]=tuple(F(int(v),10**15) for v in x[1:3])
    assert set(rows)=={'proxy_G','V2','sine2','bridge2','N_proxy','projection_error','T_proxy',*(f'third_corr{i}' for i in range(7))}
    assert all(l<=h for l,h in rows.values())
    n=rows['N_proxy']; pn=tuple(map(F,second[k]['proxy_N']))
    assert 0<=n[0]<=n[1]<F(11,20)**2 and max(n[0],pn[0])<=min(n[1],pn[1])
    eta=rows['projection_error'][1];assert rows['projection_error'][0]<=0<=eta<F(1,10)
    t=rows['T_proxy']
    c=[[-2*w*v for w,v in product(weights[i],rows[f'third_corr{i}'])] for i in range(7)]
    recomputed=(sum(map(min,c)),sum(map(max,c)))
    assert max(t[0],recomputed[0])<=min(t[1],recomputed[1])
    eps=F(old[k]['proxy_error_upper'])
    delta=s*eps+F(1,10**13)+eta
    norm=F(11,20)+eta;assert norm<1
    transport=s*(2*norm*delta+delta**2)+F(1,10**12)
    T=(t[0]-transport,t[1]+transport)
    G,M,K=(tuple(map(F,old[k][x])) for x in ('true_G','true_M','K'))
    N=tuple(map(F,second[k]['true_N']))
    a,alpha,beta=(F(targets[k][x]) for x in ('cubic_a','variational_alpha','variational_beta'))
    def upper(gg,mm,nn,tt):
        h=s*s*gg-nn
        return (g*gg-mm+(g*h+s*s*mm-tt)/a**2-2*h/a)/(ell*u)
    def lower(gg,mm,nn,tt):
        return 2*alpha*gg+2*beta*mm-alpha*alpha*(g*gg+mm)-2*alpha*beta*(g*mm+nn)-beta*beta*(g*nn+tt)
    corners=list(product(G,M,N,T))
    hi=min(max(upper(*x) for x in corners),F(second[k]['inverse_bounds'][1]))
    lo=max(min(lower(*x) for x in corners),F(second[k]['inverse_bounds'][0]))
    assert lo<=hi
    status='PASS_WITNESS' if hi<K[0] else 'FAIL_COMPARISON_WITNESS' if lo>K[1] else 'UNDECIDED'
    result[k]={'T_proxy':list(map(str,t)),'projection_error':str(eta),'transport':str(transport),'true_T':list(map(str,T)),'inverse_bounds':list(map(str,(lo,hi))),'K':list(map(str,K)),'status':status}
    print(k,'T',*[float(x) for x in T],'inverse',float(lo),float(hi),'K',*[float(x) for x in K],status)
(logs/'third_result.json').write_text(json.dumps(result,indent=2)+'\n')
