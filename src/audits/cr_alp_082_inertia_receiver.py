#!/usr/bin/env python3
"""CR-ALP-082: targeted first moments, coordinate and full-inertia receivers.

For a_i >= (R*L^-1 R)_ii and positive d_i with sum a_i/d_i <= 1,
B=R*L^-1 R <= diag(d_i). Exact positive inertia of K_lower-diag(d_i)
gives a positive subspace, without requiring a coordinate prefix.
Finite-Gram variant uses the same argument on B-G_f/u.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_075_shifted_energy import lower_matrix,positive,S,ceil
from cr_alp_078_moment_receiver import moments,sha,candidate
from cr_alp_079_inverse_gram_floor import lower
from cr_alp_081_targeted_receiver import TARGET as PREV
TARGET={'1':[27,6,16,19,21,18],'2':[20,8,21,18,28,22]}

def inertia(A):
    """Exact symmetric congruence elimination, including zero diagonal case."""
    B=[[F(v) for v in r] for r in A];p=m=z=0
    while B:
        n=len(B);k=next((i for i in range(n) if B[i][i]),None)
        if k is None:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if B[i][j]),None)
            if pair is None:z+=n;break
            i,j=pair
            # e_i -> e_i+e_j, a nonsingular basis change.
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

def test_inertia(H,aa,attempts,best):
    import numpy as np
    n=len(H);h=np.asarray(H,dtype=float)/S
    def grid(a):return max(1,ceil(a))
    def exact(d):
        assert len(d)==n and min(d)>0
        theta=sum(a*S/F(x) for a,x in zip(aa,d));assert theta<=1
        A=[[v-(d[i] if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(H)]
        return {'inertia':inertia(A),'diagonal_grid':list(map(str,d)),'theta_upper':str(theta)}
    proposals=[]
    for r,v in attempts.items():
        if not v or not v.get('positive_pivots'):continue
        r=int(r);P=list(map(int,v['diagonal_grid']));assert len(P)==r
        T=sum(aa[i]*S/F(P[i]) for i in range(r));assert T>=0
        for eta in ([F(0)] if r==n else [F(1,1000),F(1,100),F(3,100),F(1,10),F(1,5),F(1,3),F(1,2)]):
            d=[grid(F(x)*T/(1-eta)) for x in P]
            if r<n:d.extend(grid(aa[i]*S*(n-r)/eta) for i in range(r,n))
            proposals.append(d)
    # Additional exact diagonal majorizers derived from positive rational weights.
    for power in (0,.25,.5,.75,1,1.5,2):
        w=[max(1,int((max(float(a),1e-20)**power)*1e12)) for a in aa]
        W=sum(w);proposals.append([grid(a*S*W/F(v)) for a,v in zip(aa,w)])
    scored=[]
    for d in proposals:
        eig=np.linalg.eigvalsh(h-np.diag(np.asarray(d,dtype=float)/S))
        score=int(np.count_nonzero(eig>0))
        scored.append((score,float(eig[-min(n,best+1)]),d))
    # Floating values order the proposals only. Exact congruence certifies rank.
    result=None
    for score,margin,d in sorted(scored,key=lambda t:(t[0],t[1]),reverse=True)[:12]:
        row=exact(d)
        if result is None or row['inertia'][0]>result['inertia'][0]:result=row
    return result

def main(p,bank,oldpath,r70path,fullpath,prevpath,oldlogs,prevlogs,newlogs,out):
    old=json.loads(oldpath.read_text());assert old['rank']==69
    r70=json.loads(r70path.read_text());assert r70['status']=='PASS_CONDITIONAL_RANK70' and F(r70['theta_upper'])<1
    full=json.loads(fullpath.read_text());assert full['status']=='CONDITIONAL_FULL85_DIAGNOSTIC'
    prev=json.loads(prevpath.read_text());assert prev['status']=='CONDITIONAL_TARGETED_MOMENT_RECEIVER' and prev['rank']==72
    E={(e['parity'],e['coordinate_zero_based']):e for e in json.loads((bank/'source_bank.json').read_text())['entries']}
    assert len(E)==85
    for e in E.values():assert sha(bank/e['file'])==e['sha256'] and F(e['norm2'])==1
    report={};total=0
    for k,h in HASHES.items():
        ss=full['sectors'][k];J=ss['coordinate_order'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        cols=old['sectors'][k]['columns'];assert [x['coordinate'] for x in cols]==J[:35]
        a={}
        for x in cols:
            j=x['coordinate'];assert x['input_sha256']==E[k,j]['sha256']
            a[j]=F(x['norm_upper'])**2/ELL
        if k=='2':
            for j,v in r70['refined_inverse_upper'].items():a[int(j)]=min(a[int(j)],F(v))
        for j in J[35:]:
            q=moments(oldlogs/f'mode_{k}_{j}.log',E[k,j]);assert q==ss['new_moments'][str(j)]
            a[j]=F(q['inverse_upper'])
        for j in PREV[k]:
            q=moments(prevlogs/f'mode_{k}_{j}.log',E[k,j]);a[j]=min(a[j],F(q['inverse_upper']))
        ps=prev['sectors'][k];assert ps['best_conditional_rank']==36
        assert {str(j):str(v) for j,v in a.items()}==ps['inverse_bounds']
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        # Independently verify prior direct rank 36 diagonal against reconstructed bounds.
        P=list(map(int,ps['attempts']['36']['direct']['diagonal_grid']))
        H36=lower_matrix(w,J[:36]);assert len(P)==36 and min(P)>0
        assert positive([[v-(P[i] if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(H36)])
        assert sum(a[j]*S/F(d) for j,d in zip(J[:36],P))<1
        refined={}
        assert set(TARGET[k])<=set(J[:35]) and not set(TARGET[k])&set(PREV[k])
        for j in TARGET[k]:
            row=moments(newlogs/f'mode_{k}_{j}.log',E[k,j]);before=a[j];a[j]=min(before,F(row['inverse_upper']))
            refined[str(j)]={**row,'prior_inverse_upper':str(before),'retained_inverse_upper':str(a[j])}
        best=36;attempts={};direct={'36':ps['attempts']['36']['direct']};gram={}
        for r in range(n,36,-1):
            I=J[:r];v=candidate(lower_matrix(w,I),[a[j] for j in I]);direct[str(r)]=v
            H2,c=lower(w,I,a);vv=candidate(H2,c);gram[str(r)]=vv
            attempts[str(r)]={'direct':v,'inverse_gram':vv}
            if v['pass_conditional'] or vv['pass_conditional']:best=max(best,r)
        H=lower_matrix(w,J);H2,c=lower(w,J,a)
        ii=test_inertia(H,[a[j] for j in J],direct,best)
        jj=test_inertia(H2,c,gram,best)
        rank=max(best,ii['inertia'][0] if ii else 0,jj['inertia'][0] if jj else 0)
        total+=rank
        report[k]={'coordinate_order':J,'refinements':refined,'inverse_bounds':{str(j):str(v) for j,v in a.items()},
          'attempts':attempts,'best_prefix_rank':best,'full_direct_inertia':ii,'full_inverse_gram_inertia':jj,'best_conditional_rank':rank}
    result={'status':'CONDITIONAL_TARGETED_INERTIA_RECEIVER','rank':total,'conditional_index_upper':85-total,'sectors':report,
      'premises':'Authenticated original fixed-collar bank, inherited continuum and numerical enclosures, finite orthogonal projection identity for Gram variant. Full matrix inertia by exact rational congruence. No unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    assert inertia([[0,1],[1,0]])==[1,1,0]
    assert inertia([[1,1],[1,1]])==[1,0,1]
    assert inertia([[-2,0,0],[0,0,0],[0,0,3]])==[1,1,1]
    main(*map(Path,sys.argv[1:]))
