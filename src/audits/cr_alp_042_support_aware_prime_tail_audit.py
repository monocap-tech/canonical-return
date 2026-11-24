#!/usr/bin/env python3
"""Exact arithmetic for CR-ALP-042. No continuum/Lean certification.

Positive paired translations, enlarged collar (1,10) in exp coordinate.
Weights are rational upper enclosures of the original seven von Mangoldt
weights.  Step-cell recursion computes the essential row mass of T^10.
"""
from fractions import Fraction as F
from bisect import bisect_right
from math import isqrt, factorial
from pathlib import Path
import json
A=F(11,10); SCALE=10**6
POWERS={2:2,3:3,4:2,5:5,7:7,8:2,9:3}

def log12(z, count=40):
    assert 1<=z<=2
    u=(z-1)/(z+1); t=u; lo=F(0)
    for k in range(count):
        lo+=2*t/(2*k+1); t*=u*u
    return lo, lo+2*t/((2*count+1)*(1-u*u))

def log_int(n):
    k=n.bit_length()-1
    l,h=log12(F(n,2**k)); l2,h2=log12(F(2))
    return l+k*l2,h+k*h2

def ceil(x):return -((-x.numerator)//x.denominator)

assert log_int(9)[1]<2*A<log_int(10)[0]
assert sum((F(37,20)**k/factorial(k) for k in range(7)),F(0))>F(44,7)
weights={}
for n,p in POWERS.items():
    root=F(isqrt(n*10**24),10**12)
    assert 0<root and root**2<=n
    weights[n]=ceil(SCALE*log_int(p)[1]/root)
assert sum(weights.values())<SCALE*F(177,50)
print('PASS exact source weights (upper numerators / 10^6):',weights)
# Step cells (b_i,b_{i+1}); exceptional breakpoints measure zero.
# Each path contributes its product of rational upper weights, without
# merging or dropping different original prime powers.
b=[F(1),F(10)]; v=[1]; stages=[]
for k in range(1,11):
    B=sorted({F(1),F(10)}|{z for t in b for n in POWERS for z in (t*n,t/n) if 1<z<10})
    V=[]
    for l,h in zip(B,B[1:]):
        y=(l+h)/2; val=0
        for n,w in weights.items():
            for z in (n*y,y/n):
                if 1<z<10:
                    i=bisect_right(b,z)-1
                    assert b[i]<z<b[i+1]
                    val+=w*v[i]
        V.append(val)
    b,v=B,V
    stages.append({'power':k,'cells':len(v),'max_row_mass_numerator':str(max(v)),'denominator':str(SCALE**k)})
    print('power',k,'cells',len(v),'max row mass ~',float(F(max(v),SCALE**k)))
assert len(v)==2189
c=F(max(v),SCALE**10)
assert c<F(21,8)**10
print('PASS ||positive original seven-prime translation sum|| < 21/8 via tenth power')
# Two independent elementary bounds on leakage and Gamma, imported from
# the classical CR-ALP-026/027 analytic proof. pi>3, pi<22/7.
cutoffs=[]
for N,R,target in [(172,117,F(1,4)),(200,135,F(1,3)),(242,165,F(3,5)),(300,200,F(3,4)),(15000,10000,F(47,10))]:
    assert 3*(N+1)>=4*A*R and R>=2
    d=F(4096)*A**5/(45*3**7*N**5)*(F(R**5,5)+F(R**3,6)+F(R,16))
    # Round the exact log lower down to 6 decimal places; leakage up.
    ll=log_int(R)[0]; L=F(ll.numerator*10**6//ll.denominator,10**6)-F(37,20)-F(2,R)
    d=F(ceil(d*10**9),10**9)
    assert L>-6 and d<1
    mu=L*(1-d)-6*d-F(21,8)
    assert mu>target,(N,mu,target)
    cutoffs.append({'N':N,'R':R,'leakage_upper':str(d),'Gamma_high_lower':str(L),'tail_lower':str(mu),'certified_margin':str(target)})
    print(f'PASS N={N} R={R}: original signed tail > {target} (rational lower ~ {float(mu):.10f})')
# Reusing the two *pinned CR-ALP-041* energy-conditioned trace bounds.
# 2|<Tp,q>| <= 5||Q Tp||^2 + (1/5)||q||^2.
assert 5*F(13,125)<1 and 5*F(79,1000)<1
assert F(47,10)-F(1,5)==F(9,2)
print('PASS CR-ALP-041 old-cutoff mixed theorem strengthened to tail margin 9/2')
result={'stage':'CR-ALP-042','source_weight_upper_numerators':weights,'weight_scale':SCALE,'power_steps':stages,'prime_norm_upper':'21/8','cutoffs':cutoffs,'old_E14_Q15000_tail_margin':'9/2','not_claimed':'E14 plus Q172 mixed positivity; E172 finite or Schur certificate; full collar; Lean; RH'}
Path(__file__).with_name('CR-ALP-042_exact_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('CR-ALP-042 exact arithmetic PASS; continuum proofs conditional; whole-collar positivity OPEN')
