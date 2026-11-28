"""Exact endpoint replay for the two frozen physical overlap witnesses.
Requires the analytic forward-error justification in VALIDATION.md and the
saved validated.cpp output. Does not independently execute quadrature.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
prep=json.loads((p/'prepare_result.json').read_text())
ell=F(51,200); u=F(1101,200); gamma=F(72,25); prime=F(21,8)
assert ell+u==2*gamma
out={}
for k in ('1','2'):
    a=prep[k]; norm2=F(a['norm2']); assert 0<norm2<F(1000001,1000000)
    raw=(p.parent/'cr_alp_056'/f'result_extended_{k}.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==a['witness_sha256']
    w=json.loads(raw); z=w['comparison_witness']['z']; assert F(sum(x*x for x in z),10**24)==norm2
    lo=hi=0
    for i,row in enumerate(w['K']):
        for j,(l,h) in enumerate(row):
            x,y=int(l)*z[i]*z[j],int(h)*z[i]*z[j];lo+=min(x,y);hi+=max(x,y)
    K=(F(lo,int(w['scale'])*10**24), F(hi,int(w['scale'])*10**24))
    assert K==(F(a['K_lower']),F(a['K_upper']))
    rows={}
    for line in (p/f'validated_{k}.log').read_text().splitlines():
        v=line.split()
        if len(v)==5 and v[-2:] == ['/', '1e15']:
            rows[v[0]]=tuple(F(int(x),10**15) for x in v[1:3])
    assert len(rows)==10 and all(l<=h for l,h in rows.values())
    G,M,A=map(rows.get,['G','M','penalty'])
    assert 0<G[0]<=G[1]<F(4,25) and M[1]<0
    # Recompute the affine penalty from the independently displayed outward endpoints.
    AA=((gamma*G[0]-M[1])/(ell*u),(gamma*G[1]-M[0])/(ell*u))
    assert AA[0]<=A[1] and A[0]<=AA[1]
    # Exact shift-mesh perturbation allowance is 1e-13 in L2; the source
    # approximation is <5e-9 times ||z||, and the removed groups are bounded in prepare.py.
    eps=F(a['discarded_group_norm_rational_upper'])+F(500001,10**14)+F(1,10**13)
    assert F(500001,10**14)**2>F(1,200000000)**2*norm2
    # 4/10 is a rigorous norm upper bound for the integrated projected proxy.
    eG=F(4,5)*eps+eps*eps
    GT=(G[0]-eG,G[1]+eG)
    MT=(M[0]-prime*eG-F(1,10**13),M[1]+prime*eG+F(1,10**13))
    # Affine penalty on the true residual, using both moment transport bounds.
    AT=((gamma*GT[0]-MT[1])/(ell*u),(gamma*GT[1]-MT[0])/(ell*u))
    diff=(K[0]-AT[1],K[1]-AT[0]); assert diff[1]<-F(21,100)
    # The exact inverse remains unresolved: Cauchy--Schwarz lower bound.
    inv_lower=GT[0]**2/(gamma*GT[1]+MT[1]); assert inv_lower<K[0]
    out[k]={key:[str(y) for y in x] for key,x in [('K',K),('proxy_G',G),('proxy_M',M),('true_G',GT),('true_M',MT),('affine_penalty',AT),('K_minus_affine',diff)]}
    out[k]['proxy_error_upper']=str(eps);out[k]['inverse_lower']=str(inv_lower)
    # Degree-two upper polynomial, (t-ell)(t-3)^2/(9 ell t) >=0.
    # In the P_Q coordinate it is (P_Q^2 - 99/200 P_Q + 81/100 I)/(459/200).
    tol=F(459,200)*K[0]+F(99,200)*MT[0]-F(81,100)*GT[1]
    target={'1':F(138,1000),'2':F(66,1000)}[k]
    assert 0<target<tol
    out[k]['quadratic_prime_squared_moment_sufficient_threshold']=str(target)
    out[k]['quadratic_threshold_available_lower']=str(tol)
    print(k,'K',*[float(x) for x in K],'G',*[float(x) for x in GT],'M',*[float(x) for x in MT], 'affine',*[float(x) for x in AT], 'K-affine',*[float(x) for x in diff], 'inverse lower',float(inv_lower))
(p/'verification_result.json').write_text(json.dumps(out,indent=2)+'\n')
