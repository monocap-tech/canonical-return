"""Exact rational transport of the new saved second-moment quadrature.
A comparison failure is NOT a negative original Weil witness.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import sys,json
p=Path(sys.argv[1]); logs=Path(sys.argv[2]); old=json.loads((p/'cr_alp_058'/'verification_result.json').read_text())
ell,u,gamma,s=map(F,('0.255','5.505','2.88','2.625'))
result={}
for k in ('1','2'):
    rows={}
    for line in (logs/f'second_{k}.log').read_text().splitlines():
        x=line.split()
        if len(x)==5 and x[-2:]==['/','1e15']: rows[x[0]]=tuple(F(int(v),10**15) for v in x[1:3])
    assert set(rows)=={'proxy_G','V2','sine2','bridge2','N_proxy'}
    assert all(l<=h for l,h in rows.values())
    assert 0<rows['proxy_G'][0]<=rows['proxy_G'][1]<F(4,25)
    assert rows['N_proxy'][0]>=0
    # Recompute scalar orthogonal-projection identity from outward displayed endpoints.
    v,a,b,n=map(rows.get,('V2','sine2','bridge2','N_proxy'))
    recon=(max(F(0),v[0]-a[1]-b[1]),v[1]-a[0]-b[0])
    assert max(recon[0],n[0])<=min(recon[1],n[1])
    prev=old[k]; G,M,K=(tuple(map(F,prev[x])) for x in ('true_G','true_M','K'))
    eps=F(prev['proxy_error_upper'])
    # ||Q P r - Q P_rat r_proxy|| <= s eps + shift error.
    # Bound ||Q P_rat r_proxy|| via actual P and shift error; no need to
    # assert that P_rat has the same operator norm.
    shift=F(1,10**13); delta=s*eps+shift
    transport=2*(s*F(2,5)+shift)*delta+delta**2
    N=(max(F(0),n[0]-transport),n[1]+transport)
    a,b=map(F,('3.04','0.85') if k=='1' else ('2.83','1.20'))
    def vals(e,x):
        c=gamma**2-(2*x+e)*gamma+x*x+2*x*e; d=2*gamma-2*x-e
        return [(nn+d*mm+c*gg)/(e*x*x) for nn,mm,gg in product(N,M,G)]
    lower=min(vals(u,b));upper=max(vals(ell,a))
    # Also retain the affine upper and Jensen lower if tighter.
    lower=max(lower,F(prev['inverse_lower']))
    upper=min(upper,F(prev['affine_penalty'][1]))
    assert lower<=upper
    status='PASS_WITNESS' if upper<K[0] else 'FAIL_COMPARISON_WITNESS' if lower>K[1] else 'UNDECIDED'
    result[k]={'proxy_N':list(map(str,n)),'true_N':list(map(str,N)),'transport':str(transport),'inverse_bounds':[str(lower),str(upper)],'K':list(map(str,K)),'status':status}
    print(k,'N',*[float(x) for x in N],'inverse',float(lower),float(upper),'K',*[float(x) for x in K],status)
(logs/'second_result.json').write_text(json.dumps(result,indent=2)+'\n')
