#!/usr/bin/env python3
"""Fixed rational diagonal certificate; no candidate optimization."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES
from cr_alp_075_shifted_energy import lower_matrix,S
def inertia(A):
    B=[[F(v) for v in r] for r in A];p=m=z=0
    while B:
        n=len(B);k=next((i for i in range(n) if B[i][i]),None)
        if k is None:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if B[i][j]),None)
            if pair is None:z+=n;break
            i,j=pair
            for t in range(n):B[i][t]+=B[j][t]
            for t in range(n):B[t][i]+=B[t][j]
            k=i;assert B[k][k]
        if k:
            B[0],B[k]=B[k],B[0]
            for row in B:row[0],row[k]=row[k],row[0]
        d=B[0][0];p+=int(d>0);m+=int(d<0)
        B=[[B[i][j]-B[i][0]*B[0][j]/d for j in range(1,n)] for i in range(1,n)]
    assert p+m+z==len(A)
    return [p,m,z]
def main(root,prior,witness,out):
    old=json.loads(prior.read_text());assert old['status']=='CONDITIONAL_TARGETED_MOMENT_RECEIVER' and old['rank']==72
    W=json.loads(witness.read_text());assert W['parent']=='17a9486487705730dbcdac8ffb3ee708c8ecb6be'
    result={}
    for k,h in HASHES.items():
        sec=old['sectors'][k];J=sec['coordinate_order'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        s=W['sectors'][k];assert J==s['coordinates'] and int(s['scale'])==S
        a=[F(sec['inverse_bounds'][str(j)]) for j in J];assert min(a)>=0
        D=list(map(int,s['best']['diagonal']));assert len(D)==n and min(D)>0
        theta=sum(x*S/F(d) for x,d in zip(a,D));assert theta<=1
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',h)
        H=lower_matrix(w,J)
        actual=inertia([[v-(D[i] if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(H)])
        result[k]={'inertia':actual,'theta_upper':str(theta),'diagonal_grid':list(map(str,D)),
                   'matrix_scale':str(S),'coordinate_order':J,'source_matrix_sha256':h}
    rank=sum(x['inertia'][0] for x in result.values())
    ans={'status':'CONDITIONAL_FROZEN_INERTIA_RECEIVER','rank':rank,'conditional_index_upper':85-rank,'sectors':result,
         'premises':'Authenticated original fixed-collar source matrix and accepted CR-ALP-081 inverse bounds; inherited continuum and numerical enclosures. No unrestricted RH claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
    assert rank>=72
if __name__=='__main__':
    assert inertia([[0,1],[1,0]])==[1,1,0]
    assert inertia([[1,1],[1,1]])==[1,0,1]
    main(*map(Path,sys.argv[1:]))
