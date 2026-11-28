#!/usr/bin/env python3
"""CR-ALP-076 finite projected Gram plus full physical norm tail receiver.

Only valid when the pinned finite Gram is the Gram of an orthogonal
projection of the SAME physical residual as the authenticated norm logs.
Then T=R*R-Gfinite >=0; T_ii <= b_i^2-Gfinite_ii_lower.
Use weighted Cauchy to majorize T by a diagonal, and exact integer
shifted pivots for the resulting lower energy matrix.
No use of a finite lower Gram as an infinite upper Gram.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_076_finite_obstruction_inertia import interval
from cr_alp_075_shifted_energy import positive
S=10**18
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def root_up(x):
    assert x>=0
    y=x*S*S
    return F(isqrt(y.numerator//y.denominator)+1,S)
def lower_bank(w,J,include_in):
    K,G,I=w['K'],w['Gout'],w['Gin'];D=int(w['scale']);N=len(J)
    mat=[];reserve=[];diag=[]
    for i in J:
        row=[];err=0
        for j in J:
            kl,kh=interval(K,i,j);gl,gh=interval(G,i,j)
            if include_in:
                il,ih=interval(I,i,j);gl+=il;gh+=ih
            if i==j:diag.append(F(gl,D))
            l=floor((F(kl,D)-F(gh,D)/ELL)*S)
            h=ceil((F(kh,D)-F(gl,D)/ELL)*S)
            m=(l+h)//2;row.append(m);err+=max(m-l,h-m)
        mat.append(row);reserve.append(err)
    assert all(mat[i][j]==mat[j][i] for i in range(N) for j in range(N))
    for i in range(N):mat[i][i]-=reserve[i]
    return mat,diag
def compare(H,d,weights):
    assert len(d)==len(weights)==len(H) and all(x>0 for x in weights)
    assert all(t>=0 for t in d)
    # (sum sqrt(d_i)|x_i|)^2 <= (sum d_i/w_i)(sum w_i x_i^2).
    c=sum(t/w for t,w in zip(d,weights))
    M=[r[:] for r in H]
    for i,w in enumerate(weights):
        M[i][i]-=ceil(c*w*S/ELL)
    # A positive result at zero is required; a negative lower margin
    # merely means this sufficient comparison did not pass.
    l,h=-10000,10001
    def test(q):
        A=[r[:] for r in M]
        for i in range(len(A)):A[i][i]-=q*S//10000
        return positive(A)
    assert test(l)
    while h-l>1:
        q=(l+h)//2
        if test(q):l=q
        else:h=q
    return F(l,10000)
def run(p,report,out):
    data=json.loads(report.read_text())
    assert data['status']=='CONDITIONAL_PREFIX_RECEIVER'
    result={}
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        rows=data['sectors'][k]['columns']
        J=[r['coordinate'] for r in rows]
        assert len(J)==len(set(J))==35
        b=[F(r['norm_upper']) for r in rows];assert all(t>=0 for t in b)
        result[k]={}
        for include_in in (False,True):
            H,diag=lower_bank(w,J,include_in)
            d=[q*q-g for q,g in zip(b,diag)]
            assert all(t>=0 for t in d),'Physical norm contradicts finite Gram lower diagonal'
            # All choices below are exact positive rational weights.
            # The floor prevents zero weights even for zero residuals.
            a=[root_up(t) for t in d]
            choices={'uniform':[F(1)]*len(J),
                     'sqrt_tail':a,
                     'tail':[max(t,F(1,S)) for t in d],
                     'fourth_root_tail':[root_up(t) for t in a],
                     'inverse_sqrt_tail':[1/t for t in a]}
            cases={name:str(compare(H,d,v)) for name,v in choices.items()}
            best=max(map(F,cases.values()))
            result[k]['Gin_plus_Gout' if include_in else 'Gout_only']={
                'tail_trace_upper':str(sum(d)),
                'margin_lower_by_weight':cases,
                'best_margin_lower':str(best),
                'pass_conditional':best>0}
    ok=all(any(x['pass_conditional'] for x in r.values()) for r in result.values())
    ans={'status':'PASS_CONDITIONAL_RANK70' if ok else 'UNDECIDED',
         'sectors':result,'conditional_index_if_received':15 if ok else 16,
         'premises':'Authenticated 075 source, norm evidence, and orthogonal projected finite Gram identity, inherited continuum and numerical enclosures. No original negative witness or RH claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))
if __name__=='__main__':run(*(Path(x) for x in sys.argv[1:]))
