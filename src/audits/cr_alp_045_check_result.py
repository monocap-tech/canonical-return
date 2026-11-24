#!/usr/bin/env python3
"""Independent exact saved-interval determinant and acceptance check."""
from fractions import Fraction as F
from pathlib import Path
import json
from cr_alp_045_layered_gamma_audit import lower,weighted
ROOT=Path(__file__).resolve().parent
S=10**36
r=json.loads((ROOT/'CR-ALP-045_fitted_extension_result.json').read_text())
def det(a):
    a=[row[:] for row in a];prev=1
    for k in range(len(a)-1):
        pivot=a[k][k];assert pivot>0
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                v=a[i][j]*pivot-a[i][k]*a[k][j]
                assert v%prev==0;a[i][j]=v//prev
        prev=pivot
    return a[-1][-1]
for b in r:
    A=[[(int(x[0]),int(x[1])) for x in row] for row in b['gram_intervals']]
    assert A==list(map(list,zip(*A)))
    R=max(sum(h-l for l,h in row) for row in A)
    band=F(b['band_upper']);v=2*S*band;assert v.denominator==1
    C=[[int(v)*(i==j)-sum(A[i][j])-R*(i==j) for j in range(7)] for i in range(7)]
    for n in range(1,8):assert det([row[:n] for row in C[:n]])>0
    assert F(b['rho_upper'])==band+F(b['tail_upper'])
    mu,_=lower(b['N']);assert F(b['mu_lower'])==mu
    assert b['pass_sufficient']==(F(b['rho_upper'])<mu)
for N,eta,mu,margin in [(112,F(107,200),F(537,1000),F(1,500)),(123,F(1,2),F(3,5),F(1,10)),(172,F(23,100),F(9,10),F(67,100))]:
    b=[x for x in r if x['N']==N];assert len(b)==2 and all(F(x['rho_upper'])<eta<mu<F(x['mu_lower']) for x in b)
    assert mu-eta==margin
    print('PASS E14 + Q%d: exact determinants and mixed tail margin'%N,margin)
assert lower(66)[0]>F(1,100000) and lower(85)[0]>F(1,4)
top,cost=weighted(123,171);assert top>F(13,20) and cost<F(1,20)
print('PASS all saved Gram comparisons, Q66/Q85, simultaneous Q123 weights; no full collar or Lean claim')
