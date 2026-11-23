#!/usr/bin/env python3
"""CR-ALP-009 exact Gaussian-rational antiunitary phase no-go.

A model of generic finite InputS/InputSG, NOT an actual arithmetic D7
three-word theta lift. No Lean compilation, signed Weil transfer, or RH claim.
All assertions use Fraction; no floating point.
"""
from fractions import Fraction as Q
from itertools import permutations

ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(0))

def z(x=0, y=0):
    return Q(x), Q(y)

def add(x, y):
    return x[0]+y[0], x[1]+y[1]

def sub(x, y):
    return x[0]-y[0], x[1]-y[1]

def mul(x, y):
    return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]

def conj(x):
    return x[0], -x[1]

def abs2(x):
    return x[0]**2+x[1]**2

def sumz(values):
    out=ZERO
    for v in values:
        out=add(out,v)
    return out

def mm(A,B):
    return [[sumz(mul(A[i][k],B[k][j]) for k in range(3))
             for j in range(3)] for i in range(3)]

def trace(A):
    return sumz(A[i][i] for i in range(3))

def determinant(A):
    return add(
        sub(mul(A[0][0],sub(mul(A[1][1],A[2][2]),mul(A[1][2],A[2][1]))),
            mul(A[0][1],sub(mul(A[1][0],A[2][2]),mul(A[1][2],A[2][0])))),
        mul(A[0][2],sub(mul(A[1][0],A[2][1]),mul(A[1][1],A[2][0]))))

a,b,c = Q(1,5),Q(1,6),Q(1,7)
S=a*a+b*b+c*c
assert S == Q(3889,44100) and 1-S == Q(40211,44100)

def gram(sign):
    return [[ONE,z(a),z(b)],
            [z(a),ONE,z(0,sign*c)],
            [z(b),z(0,-sign*c),ONE]]

Gplus,Gminus=gram(1),gram(-1)
for sign,G in [(1,Gplus),(-1,Gminus)]:
    assert all(G[j][i]==conj(G[i][j]) for i in range(3) for j in range(3))
    # Sylvester positivity; rank=3.
    assert 1-a*a>0 and determinant(G)==z(1-S) and 1-S>0
    # Literal generic InputS/InputSG: weights 1, D=D*=3, all defects 0.
    assert all(G[i][i]==ONE for i in range(3))
    assert max(abs2(G[i][j]) for i in range(3) for j in range(3)
               if i!=j)<=Q(1,3)
    assert trace(G)==z(3)
    assert trace(mm(G,G))==z(3+2*S)
    assert trace(mm(mm(G,G),G))==z(3+6*S)
    assert mul(mul(G[0][1],G[1][2]),G[2][0])==z(0,sign*a*b*c)
    # The exact characteristic polynomial is (lambda-1)^3-S(lambda-1).
    for lam in [-2,0,1,2,5]:
        A=[[sub(z(lam) if i==j else ZERO,G[i][j])
            for j in range(3)] for i in range(3)]
        assert determinant(A)==z((Q(lam)-1)**3-S*(Q(lam)-1))
    # All 8 independent real sign gauges leave cyclic orientation fixed.
    for bits in range(8):
        eps=[z(1 if not (bits & (1<<j)) else -1) for j in range(3)]
        h=[[mul(mul(conj(eps[i]),G[i][j]),eps[j])
            for j in range(3)] for i in range(3)]
        assert mul(mul(h[0][1],h[1][2]),h[2][0])==z(0,sign*a*b*c)

# Equal complete spectral data under complex conjugation, not only 2 moments.
assert Gminus==[[conj(v) for v in row] for row in Gplus]
P,N=Gplus,Gminus
for _ in range(1,15):
    assert trace(P)==trace(N) and trace(P)[1]==0
    P,N=mm(P,Gplus),mm(N,Gminus)

# Unequal labelled edge magnitudes leave only identity label permutation.
for perm in permutations(range(3)):
    if all(abs2(Gplus[i][j])==abs2(Gminus[perm[i]][perm[j]])
           for i in range(3) for j in range(i+1,3)):
        assert perm==(0,1,2)

print('CR-ALP-009 GAUSSIAN-RATIONAL PASS')
print('S =',S,' determinant =',1-S)
print('Spectrum both: 1, 1±sqrt(S); all spectral moments coincide')
print('Two labelled cyclic traces: +i/210 and -i/210')
print('Generic InputS/InputSG tests pass; actual arithmetic D7 remains uninstantiated')
