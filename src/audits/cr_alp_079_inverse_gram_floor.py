#!/usr/bin/env python3
"""CR-ALP-079: finite projected Gram floor for the inverse-energy matrix.

Let B=R* L^-1 R with ell I <= L <= u I, u=1101/200.
For an orthogonal finite projection P of the SAME residual,
B >= G_f/u with G_f=R* P R. Consequently B-G_f/u is PSD.
Use coordinate inverse upper bounds a_i and finite Gram diagonal lower
endpoints to bound its diagonal. Retain the finite Gram exactly in the
energy comparison. No promotion from finite Gram alone.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_075_shifted_energy import S,floor,ceil
from cr_alp_078_moment_receiver import candidate
U=F(1101,200)
def interval(A,i,j):
    a,b=map(int,A[i][j]);c,d=map(int,A[j][i])
    l,h=max(a,c),min(b,d);assert l<=h
    return l,h
def lower(w,J,a):
    K,G=w['K'],w['Gout'];D=int(w['scale'])
    H=[];res=[];c=[]
    for i in J:
        row=[];err=0
        for j in J:
            kl,kh=interval(K,i,j);gl,gh=interval(G,i,j)
            l=floor((F(kl,D)-F(gh,D)/U)*S)
            h=ceil((F(kh,D)-F(gl,D)/U)*S)
            assert l<=h
            m=(l+h)//2;row.append(m);err+=max(m-l,h-m)
            if i==j:
                # Nonnegativity of the exact finite Gram diagonal is known.
                ci=a[i]-max(F(0),F(gl,D))/U
                assert ci>=0,('inverse bound inconsistent with finite floor',i)
                c.append(ci)
        H.append(row);res.append(err)
    assert all(H[i][j]==H[j][i] for i in range(len(J)) for j in range(len(J)))
    for i,e in enumerate(res):H[i][i]-=e
    return H,c
def main(p,oldpath,rank70path,fullpath,out):
    old=json.loads(oldpath.read_text());assert old['rank']==69
    r70=json.loads(rank70path.read_text());assert r70['status']=='PASS_CONDITIONAL_RANK70'
    assert F(r70['theta_upper'])<1
    full=json.loads(fullpath.read_text());assert full['status']=='CONDITIONAL_FULL85_DIAGNOSTIC'
    report={};total=0
    for k,h in HASHES.items():
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        ss=full['sectors'][k];J=ss['coordinate_order'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        cols=old['sectors'][k]['columns'];assert [x['coordinate'] for x in cols]==J[:35]
        a={x['coordinate']:F(x['norm_upper'])**2/ELL for x in cols}
        if k=='2':
            for j,v in r70['refined_inverse_upper'].items():a[int(j)]=min(a[int(j)],F(v))
        assert sorted(map(int,ss['new_moments']))==sorted(J[35:])
        for j,v in ss['new_moments'].items():a[int(j)]=F(v['inverse_upper'])
        assert all(v>=0 for v in a.values())
        best=int(ss['best_conditional_rank']);assert 35<=best<=n
        attempts={}
        for r in range(n,best,-1):
            H,c=lower(w,J[:r],a)
            v=candidate(H,c)
            attempts[str(r)]=v
            if v['pass_conditional']:best=r;break
        total+=best
        report[k]={'best_conditional_rank':best,'attempts':attempts}
    result={'status':'CONDITIONAL_INVERSE_GRAM_FLOOR',
            'rank':total,'conditional_index_upper':85-total,'sectors':report,
            'premises':'Authenticated 056/057/075/077/078 evidence; same-residual finite orthogonal projection identity, original L spectral bounds and inherited continuum/numerical enclosures. No unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
