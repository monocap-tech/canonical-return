#!/usr/bin/env python3
"""Independent exact check of saved fitted-band matrix intervals and acceptance ledger.
Does not reconstruct the coefficient sweep or certify the continuum tail identities.
"""
from fractions import Fraction as F
from pathlib import Path
import json
S=10**36
ROOT=Path(__file__).resolve().parent
D=json.loads((ROOT/'CR-ALP-044_exact_result.json').read_text())
assert D['stage']=='CR-ALP-044' and F(D['sigma'])==F(1,10**14)
assert [b['parity'] for b in D['blocks']]==[1,2]
eta=F(D['eta']);mu=F(1,4);assert eta==F(23,100) and mu-eta==F(D['tail_margin'])==F(1,50)
def det(A):
    A=[r[:] for r in A];n=len(A);prev=1
    for k in range(n-1):
        pivot=A[k][k];assert pivot>0
        for i in range(k+1,n):
            for j in range(k+1,n):
                v=A[i][j]*pivot-A[i][k]*A[k][j]
                assert v%prev==0;A[i][j]=v//prev
        prev=pivot
    return A[-1][-1]
for b in D['blocks']:
    assert (b['N'],b['K'],b['rows'])==(172,30000,14914)
    band=F(b['band_operator_upper']);assert band=={1:F(9,50),2:F(4,25)}[b['parity']]
    g=[[(int(l),int(h)) for l,h in r] for r in b['fitted_band_gram_intervals']]
    assert len(g)==7 and all(len(r)==7 for r in g)
    assert all(l<=h for r in g for l,h in r)
    assert all(g[i][j]==g[j][i] for i in range(7) for j in range(7))
    rad=max(sum(h-l for l,h in r) for r in g)
    assert F(rad,2*S)<=F(b['band_radius_upper'])
    # All entries are integers in scale 2*S. Every admissible real
    # symmetric Gram lies below its midpoint plus rad/(2*S) I.
    assert (2*S*band).denominator==1
    A=[[int(2*S*band)*(i==j)-sum(g[i][j])-rad*(i==j) for j in range(7)] for i in range(7)]
    assert all(det([r[:n] for r in A[:n]])>0 for n in range(1,8))
    # Use separately upward-rounded component norm budgets. This is a
    # slightly looser independent final inequality, not a new tail proof.
    tail=(sum(F(b[k]) for k in ['log_norm_upper','prime_norm_upper','error_norm_upper']))**2
    assert band+tail<eta<mu
    assert F(b['rho_upper'])==band+F(b['tail_operator_upper'])<eta
    print('PASS parity',b['parity'],': seven exact Bareiss determinants; rounded component norm budgets below eta')
assert 172-14==158
print('PASS saved interval-matrix and rational acceptance-ledger check; conditional analytical result only')
