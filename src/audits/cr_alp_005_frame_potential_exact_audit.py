#!/usr/bin/env python3
"""CR-ALP-005: exact finite frame ambiguity and trace-rank regression.

Exact symbolic algebra ONLY, not a proof of any arithmetic theta-lift
inequality, Weil global positivity, operator-domain transfer, or RH.
"""
import sympy as sp

b = sp.symbols("b", real=True, nonnegative=True)
a = sp.Integer(9)
w = sp.Integer(1)
L = sp.log(9)
mass = a-w

# f(u)=exp(u/2) 1_[0,log9]; f_b(u)=f(u-b).
# For 0 <= b <= L the cross physical pairing is:
def overlap(shift):
    return a * sp.exp(-shift/2) - w * sp.exp(shift/2)

assert sp.simplify(overlap(sp.Integer(0))-mass)==0
assert sp.simplify(overlap(sp.log(3))-2*sp.sqrt(3))==0
assert sp.simplify(overlap(L))==0

def gram(c):
    return sp.Matrix([[mass,c],[c,mass]])

def isotropy_excess(M):
    # D=2, choose a 2-dimensional target containing the source span.
    # For b=0, the target has one additional unused direction.
    return sp.factor(2*sp.trace(M*M)/sp.trace(M)**2-1)

M_equal=gram(mass)
M_overlap=gram(2*sp.sqrt(3))
M_orth=gram(sp.Integer(0))

assert M_equal.rank()==1
assert M_overlap.rank()==2
assert M_orth.rank()==2
assert M_overlap.det()==52
assert isotropy_excess(M_equal)==1
assert isotropy_excess(M_overlap)==sp.Rational(3,16)
assert isotropy_excess(M_orth)==0

# The individual arithmetic table rows are independent of translation.
for n in [2,3,4,9,27,100]:
    row=(sp.Rational(9,n)-1) if n<=9 else sp.Integer(0)
    assert row==max(sp.Rational(9,n)-1,sp.Integer(0))

# Frame isotropy fixes neither an independently specified Hermitian
# numerator nor its sign.
B=sp.eye(2)
Q=B*B.T
H_positive=sp.eye(2)
H_indefinite=sp.diag(-1,1)
assert isotropy_excess(B.T*B)==0
assert H_positive.det()==1 and H_indefinite.det()==-1
assert Q==sp.eye(2)

# rank(Q)*tr(Q^2) >= tr(Q)^2 for PSD Q, sample exact cases.
for diagonal in [(1,1),(1,2),(0,2),(3,7),(1,0),(0,0)]:
    Q=sp.diag(*diagonal)
    if sp.trace(Q)!=0:
        assert Q.rank()*sp.trace(Q*Q)>=sp.trace(Q)**2

print("CR-ALP-005 EXACT SYMBOLIC PASS")
print("Translated-source Gram ranks: 1, 2, 2")
print("Isotropy excess for D=2: 1, 3/16, 0")
print("Scope: finite algebra only; no arithmetic (M), Lean, or RH result")
