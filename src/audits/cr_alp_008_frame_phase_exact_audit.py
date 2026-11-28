#!/usr/bin/env python3
"""CR-ALP-008 pure-rational test of phase blindness in the finite InputS schema.

Exact finite linear algebra only; this does NOT instantiate the arithmetic D7
word operators, compile Lean, prove (M), or transport a signed Weil operator.
"""
from fractions import Fraction as Q
from itertools import product


def mat(sign, r):
    return [[Q(1), r, r], [r, Q(1), sign * r], [r, sign * r, Q(1)]]


def multiply(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(3)), Q(0))
             for j in range(3)] for i in range(3)]


def trace(a):
    return sum((a[i][i] for i in range(3)), Q(0))


def det3(a):
    return (a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
           -a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
           +a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]))

for r in [Q(1,4), Q(1,8), Q(1,100), Q(1,10000)]:
    for sign in (-1, 1):
        g = mat(sign, r)
        g2 = multiply(g,g)
        g3 = multiply(g2,g)
        assert det3(g) == 1-3*r*r+2*sign*r*r*r > 0
        assert 1-r*r > 0     # Sylvester leading two-by-two minor
        assert trace(g) == 3
        assert trace(g2) == 3+6*r*r
        assert trace(g3) == 3+18*r*r+6*sign*r*r*r
        assert 3*trace(g2)/(trace(g)**2)-1 == 2*r*r
        assert g[0][1]*g[1][2]*g[2][0] == sign*r**3
        assert all(g[i][i] == 1 for i in range(3))
        assert all(g[i][j]**2 == r*r for i in range(3)
                   for j in range(3) if i != j)
        # Abstract Lean InputS finite-frame schema:
        # D=Dstar=3, ga_i=1, d_i=1, eta0=eta1=e0=0.
        assert r*r <= Q(1,3)   # off-diagonal HS squared bound
        assert 3 <= 3          # target dimension bound
        for eps in product((-1,1), repeat=3):
            h=[[eps[i]*g[i][j]*eps[j] for j in range(3)]
               for i in range(3)]
            assert h[0][1]*h[1][2]*h[2][0] == sign*r**3

r=Q(1,4)
plus=mat(1,r); minus=mat(-1,r)
assert det3(plus)==Q(27,32)
assert det3(minus)==Q(25,32)
assert trace(multiply(plus,plus))==Q(27,8)
assert trace(multiply(minus,minus))==Q(27,8)
assert trace(multiply(multiply(plus,plus),plus))==Q(135,32)
assert trace(multiply(multiply(minus,minus),minus))==Q(129,32)
assert trace(multiply(multiply(plus,plus),plus))-trace(multiply(multiply(minus,minus),minus))==Q(3,16)
print('CR-ALP-008 PURE-RATIONAL FINITE FRAME AUDIT PASS')
print('Same rank=3, diagonals, |Kij|^2=1/16, trQ=3, trQ²=27/8, Phi3=1/8')
print('Opposite cyclic trace: +1/64 versus -1/64; trQ³: 135/32 versus 129/32')
print('All 8 independent real sign gauges preserve the cyclic sign')
print('The test covers the abstract InputS shape, not the actual theta-lift D7 family')
