#!/usr/bin/env python3
"""CR-ALP-073: deterministic rank-20 source selection, no new physical certificate.

Input is the hash-verified extracted 055-058 archive. Retain 070's exact
first source; select 19 even coordinate directions using only K_L and
upper bounds on their correction energy. Output an exact robust block
receiving budget, not a claim about unknown physical residual norms.
"""
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
D=10**12
S=10**30
def up(x):
    assert x>=0;y=x*S*S
    return F(isqrt(y.numerator//y.denominator)+1,S)
def read(p,h):
    raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()==h,p
    return json.loads(raw)
def absmax(q):return max(map(abs,q))
def run(package,out):
    w=read(package/'cr_alp_056/result_extended_2.json','4017c7a656e43fe9fa97da91c37f35a9eda5a78dd3efaa2863b40de71177ac36')
    trial=read(package/'cr_alp_057/trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    z=list(map(int,w['comparison_witness']['z']));n=len(z);assert n==42
    ks=int(w['scale']);K=w['K'];assert len(K)==n and all(len(row)==n for row in K)
    def entry(i,j):
        a,b=(F(int(x),ks) for x in K[i][j])
        c,d=(F(int(x),ks) for x in K[j][i])
        l,h=max(a,c),min(b,d);assert l<=h
        return l,h
    def contraction(v):
        l=h=F(0)
        for i,x in enumerate(z):
            a,b=entry(i,v);u=F(x,D)
            l+=min(a*u,b*u);h+=max(a*u,b*u)
        return l,h
    cross=[absmax(contraction(i)) for i in range(n)]
    ts=int(trial['scale']);assert ts==10**100
    t=trial['sectors']['2'];assert t['ids']==list(range(2,851,2))
    assert all(len(row)==n for row in t['D'])
    A=F(555992307,125000000)-F(72,25);assert A>0
    ai=[A*sum(absmax(tuple(F(int(x),ts) for x in row[i]))**2 for j,row in zip(t['ids'],t['D']) if j>85) for i in range(n)]
    # Anchor correction bound uses the actual contracted finite precursor.
    a1=F(0)
    for j,row in zip(t['ids'],t['D']):
        if j<=85:continue
        l=h=0
        for (a,b),x in zip(row,z):
            a,b=int(a),int(b);assert a<=b
            l+=min(a*x,b*x);h+=max(a*x,b*x)
        a1+=max(abs(F(l,ts*D)),abs(F(h,ts*D)))**2
    a1*=A
    m=F(61214171350249519,10**18);B1=F(16059528681514303,62500000000000000);dmin=F(26535863,100000000)
    def bound(J,q):
        # K_L easy block >= lambda I via symmetric Gershgorin.
        lam=min(entry(i,i)[0]-sum(absmax(entry(i,j)) for j in J if j!=i) for i in J)
        a=sum(ai[i] for i in J);h=up(sum(cross[i]**2 for i in J))
        # q is a hypothetical common actual old-residual norm ceiling.
        B=sum((q*up(1/dmin)+up(ai[i]))**2 for i in J)
        c=h+up(a1*a)+up(B1*B)
        mu=lam-B;gap=m*mu-c*c
        return gap,mu,lam,B,c,a,h
    # Include the already attached second direction (index 37).
    J=[37]
    for _ in range(18):
        candidates=[]
        for i in range(n):
            if i in J:continue
            jj=J+[i]
            if not any(z[k] for k in range(n) if k not in jj):continue
            candidates.append((bound(jj,F(0))[0],-i,i))
        assert candidates
        J.append(max(candidates)[2])
    assert len(J)==19 and len(set(J))==19
    outside=[i for i in range(n) if i not in J and z[i]!=0];assert outside
    # Monotone search of a sufficient common norm threshold at 1e-5 grid.
    zero=bound(J,F(0));l=0;h=20001
    if zero[0]>0 and zero[1]>0:
        while h-l>1:
            mid=(h+l)//2;b=bound(J,F(mid,10**5))
            if b[0]>0 and b[1]>0:l=mid
            else:h=mid
    q=F(l,10**5);gap,mu,lam,B,c,a,kn=bound(J,q)
    result={'status':'CANDIDATE_ONLY','rank':20,'anchor_numerators':z,'source_denominator':D,
      'coordinate_indices_zero_based':J,'independence_witness_outside_index':outside[0],
      'coordinate_correction_energy_uppers':{str(i):str(ai[i]) for i in J},
      'anchor_correction_energy_upper':str(a1),'easy_K_lower':str(lam),
      'K_cross_norm_upper':str(kn),'sum_correction_energy_upper':str(a),
      'common_actual_residual_norm_threshold':str(q),'sum_inverse_upper_at_threshold':str(B),
      'easy_S_lower_at_threshold':str(mu),'cross_norm_upper_at_threshold':str(c),
      'block_test_margin_at_threshold':str(gap),
      'conditional_arithmetic_pass':mu>0 and gap>0,
      'required_new_evidence':'Each selected actual residual norm <= threshold, with inherited analytical premises. No index promotion.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))
