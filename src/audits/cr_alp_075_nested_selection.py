#!/usr/bin/env python3
"""CR-ALP-075 nested expansion. Exact source selection only.
Reuses original SHA-pinned 056 matrices and the existing rank-20 prefix.
The 30+30 candidate is not a physical positivity certificate.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
PREFIX={'1':[0,34,25,35,1,22,32,39,2,40],
        '2':[37,0,25,1,38,2,27,3,29,4]}
def select(w,k):
    K=w['K'];S=int(w['scale']);n=len(K)
    assert n==(43 if k=='1' else 42)
    M=[]
    for i in range(n):
        row=[]
        for j in range(n):
            a,b=map(int,K[i][j]);c,d=map(int,K[j][i])
            a,b=max(a,c),min(b,d);assert a<=b
            row.append((a,b))
        M.append(row)
    def lower(J):
        return F(min(M[i][i][0]-sum(max(map(abs,M[i][j])) for j in J if i!=j) for i in J),S)
    J=[];levels={}
    for r in range(1,n+1):
        J.append(max((lower(J+[j]),-j,j) for j in range(n) if j not in J)[2])
        if r>=10:
            lam=lower(J);levels[str(r)]={'K_lower':str(lam),'residual_sum_budget':str(ELL*lam)}
    assert J[:10]==PREFIX[k]
    lam=lower(J[:30]);assert lam>0
    q2=ELL*lam/60;Q=10**12;t=q2*Q*Q
    q=F(isqrt(t.numerator//t.denominator),Q)
    return {'rank':30,'indices_zero_based':J[:30],
            'old_prefix':J[:10],'new_coordinates':J[10:30],
            'complete_order':J,'prefix_levels':levels,
            'K_lower':str(lam),'common_actual_residual_norm_ceiling':str(q),
            'conditional_gate_usable':q>0}
def main(p,out):
    result={'status':'SOURCE_SELECTION_ONLY','rank':60,
        'coercivity_lower':str(ELL),'conditional_gate_usable':True,
        'sectors':{k:select(checked(p/'cr_alp_056'/f'result_extended_{k}.json',h),k) for k,h in HASHES.items()},
        'conditional_index_target':25}
    assert all(v['conditional_gate_usable'] for v in result['sectors'].values())
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
