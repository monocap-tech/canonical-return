#!/usr/bin/env python3
"""Exact rational arithmetic for independent layered Gamma lower bounds.
Continuum identities inherited from CR-ALP-026/027/042; no Lean claim.
"""
from fractions import Fraction as F
from math import factorial
import json
from pathlib import Path
A=F(11,10); PI=F(314159,100000)
# Machin identity: pi = 16 atan(1/5)-4 atan(1/239), alternating bounds.
def atan(x,n):
    return sum(((-1)**k*x**(2*k+1)/ (2*k+1) for k in range(n)),F(0))
assert 16*atan(F(1,5),8)-4*atan(F(1,239),7)>PI
assert 16*atan(F(1,5),7)-4*atan(F(1,239),8)<F(22,7)
assert sum((F(37,20)**k/factorial(k) for k in range(7)),F(0))>F(44,7)
def log12(x,n=24):
    u=(x-1)/(x+1)
    return 2*sum((u**(2*k+1)/(2*k+1) for k in range(n)),F(0))
LOG2=log12(F(2))
def loglo(n):
    k=n.bit_length()-1; v=k*LOG2+log12(F(n,2**k));return F(v*10**9//1,10**9)
def ceil(x):return -((-x.numerator)//x.denominator)
def ant(r):return F(r**5,5)+F(r**3,6)+F(r,16)
def lower(N):
    J=N+1; nu=J*PI/(2*A)
    # each parity tail sum <= J^-6 + 1/(10J^5).
    c=4/(A*PI)*(2*A/PI)**6*(F(1,J**6)+F(1,10*J**5))
    d=F(0); low=F(-6); prev=F(-6); r0=0
    # each incremental low-band mass bounded by the exact integrated
    # numerator times the worst denominator on its unit interval.
    for r in range(1,int(nu)):
        inc=c*(ant(r)-ant(r-1))/(1-(r/nu)**2)**2
        d+=F(ceil(inc*10**12),10**12)
        if r>=2:
            L=loglo(r)-F(37,20)-F(2,r)
            assert L>prev
            low+=(L-prev)*max(F(0),1-d)
            prev=L
        if d>=1:
            r0=r;break
    return low-F(21,8),r0
def weighted(N,T,s=2):
    J=N+1; nu=J*PI/(2*A)
    assert 0<=s<5 and 2<=T<nu
    c=4/(A*PI)*(2*A/PI)**6*(F(1,J**6)+F(1,2*(5-s)*J**5))
    d=F(0); cost=F(0); prev=F(-6)
    for r in range(1,T+1):
        d+=F(ceil(c*(ant(r)-ant(r-1))/(1-(r/nu)**2)**2*10**12),10**12)
        if r>=2:
            L=loglo(r)-F(37,20)-F(2,r)
            assert L>prev
            cost+=(L-prev)*d;prev=L
    return prev-F(21,8),cost

if __name__=='__main__':
    # strict lower targets, chosen independently of Return.
    result={}
    for N,target in [(66,F(1,100000)),(67,F(1,100)),(85,F(1,4)),(90,F(3,10)),(100,F(2,5)),(123,F(3,5)),(150,F(4,5)),(172,F(9,10)),(242,F(6,5)),(300,F(3,2))]:
        v,r=lower(N);print(N,r,float(v),target,flush=True)
        assert v>target
        result[str(N)]={'lower':str(v),'last_radius':r,'target':str(target),'pass':v>target}
    top,cost=weighted(123,171)
    assert top>F(13,20) and cost<F(1,20)
    result['weighted_Q123']={'baseline':str(top),'penalty':str(cost),'certified_baseline':'13/20','certified_penalty':'1/20','weight':'13/20-(124/j)^2/20'}
    print('PASS simultaneous weighted Q123: 13/20 - (124/j)^2/20')
    Path(__file__).with_name('CR-ALP-045_layered_gamma_result.json').write_text(json.dumps(result,indent=2)+'\n')
