#!/usr/bin/env python3
"""CR-ALP-075 nested expansion. Exact source selection only.
Reuses original SHA-pinned 056 matrices and the existing rank-20 prefix.
The 35+35 candidate is not a physical positivity certificate.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
PREFIX={'1':[0,34,25,35,1,22,32,39,2,40],
        '2':[37,0,25,1,38,2,27,3,29,4]}
def select(w,k,target=30):
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
    assert 10<=target<=35
    lam=lower(J[:target]);assert lam>0
    q2=ELL*lam/(2*target);Q=10**12;t=q2*Q*Q
    q=F(isqrt(t.numerator//t.denominator),Q)
    return {'rank':target,'indices_zero_based':J[:target],
            'old_prefix':J[:10],'new_coordinates':J[10:target],
            'complete_order':J,'prefix_levels':levels,
            'K_lower':str(lam),'common_actual_residual_norm_ceiling':str(q),
            'conditional_gate_usable':q>0}
def main(p,out):
    # Import after select is defined: shifted_energy also imports select.
    from cr_alp_075_shifted_energy import lower_matrix,shift
    target=35
    sectors={}
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        s=select(w,k,target)
        J=s['indices_zero_based']
        for n in range(10,target+1):
            level=s['prefix_levels'][str(n)]
            g=F(level['K_lower']);q=shift(lower_matrix(w,J[:n]))
            level['gershgorin_K_lower']=str(g)
            level['shifted_K_lower']=str(q)
            level['K_lower']=str(max(g,q))
            level['residual_sum_budget']=str(ELL*max(g,q))
        lam=F(s['prefix_levels'][str(target)]['K_lower'])
        Q=10**12;t=ELL*lam*Q*Q/(2*target)
        s['K_lower']=str(lam)
        s['common_actual_residual_norm_ceiling']=str(F(isqrt(t.numerator//t.denominator),Q))
        s['conditional_gate_usable']=lam>0
        sectors[k]=s
    result={'status':'SOURCE_SELECTION_ONLY','rank':2*target,
        'coercivity_lower':str(ELL),'conditional_gate_usable':True,
        'sectors':sectors,'conditional_index_target':85-2*target}
    assert all(v['conditional_gate_usable'] for v in sectors.values())
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
