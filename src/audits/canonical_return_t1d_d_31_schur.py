#!/usr/bin/env python3
"""CR-INT-T1D-D exact collision/Gram/Schur audit.

Checks prime-power window census, ordered collision thresholds, the
31-prime four-cell-to-two-cell perturbation, pole orthogonality, the
selected-gap preservation inequality through the next spacing-two event,
and the finite-form Schur completion identity.

NOT a Lean certification, numerical proof of analytic Gamma bounds, or
an operator/complement positivity theorem.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp


def prime_power(n: int) -> bool:
    if n < 2:
        return False
    for p in range(2, n + 1):
        if n % p:
            continue
        if any(p % q == 0 for q in range(2, isqrt(p) + 1)):
            continue
        k = n
        while k % p == 0:
            k //= p
        if k == 1:
            return True
    return False


# For ell < log(11/9), n contributes at offset k log 3 only if
# 3**k * (9/11) < n < 3**k * (11/9).
CUTOFF = F(11, 9)
windows = []
for m in (3, 9, 27):
    windows.append([n for n in range(2, 80)
                    if prime_power(n) and F(9*m, 11) < n < F(11*m, 9)])
assert windows == [[3], [8, 9], [23, 25, 27, 29, 31, 32]], windows

thresholds = [
    ('29', F(29, 27)),
    ('25', F(27, 25)),
    ('8', F(9, 8)),
    ('31', F(31, 27)),
    ('23', F(27, 23)),
    ('32', F(32, 27)),
    ('11', F(11, 9)),
]
assert all(thresholds[i][1] < thresholds[i + 1][1] for i in range(len(thresholds) - 1))
# Half-support bound for t = log(9/8) through log(11/9):
assert F(11, 9) < F(9, 8) ** 2
assert F(5, 3) + F(1, 8) < 2

r = sp.sqrt(3)
s = 4 / r
X = sp.Matrix([[1, 0], [-s, 1], [1, -s], [0, 1]])
M = sp.simplify(X.T * X)
expected_M = sp.Matrix([[sp.Rational(22, 3), -8*r/3],
                        [-8*r/3, sp.Rational(22, 3)]])
assert sp.simplify(M - expected_M) == sp.zeros(2, 2)
L = sp.Matrix([[r**i for i in range(4)], [r**(-i) for i in range(4)]])
assert sp.simplify(L*X) == sp.zeros(2, 2)

# The new outermost-cell 31 interaction, also valid for 23 and 32.
D3 = sp.zeros(4)
D3[0, 3] = D3[3, 0] = 1
B3 = sp.simplify(X.T*D3*X)
assert B3 == sp.Matrix([[0, 1], [1, 0]])

m_plus = sp.simplify(M[0,0] + M[0,1])
m_minus = sp.simplify(M[0,0] - M[0,1])
# New Weil increment is -J * kappa_n * R_phi(delta_n) * B3.
plus_coefficient = sp.simplify(-1/m_plus)
minus_coefficient = sp.simplify(1/m_minus)
split_coefficient = sp.simplify(minus_coefficient-plus_coefficient)
assert sp.simplify(plus_coefficient + 3/(22-8*r)) == 0
assert sp.simplify(minus_coefficient - 3/(22+8*r)) == 0
assert split_coefficient == sp.Rational(33,73)

# The previous prime-8 (spacing-two) contribution has opposite gap sign.
D2 = sp.zeros(4)
for i in (0,1):
    D2[i,i+2] = D2[i+2,i] = 1
B2 = sp.simplify(X.T*D2*X)
assert B2 == sp.Matrix([[2,-8*r/3],[-8*r/3,2]])
split8 = sp.simplify(-(B2[0,0]-B2[0,1])/m_minus
                     + (B2[0,0]+B2[0,1])/m_plus)
assert sp.simplify(split8 + 64*r/73) == 0


def log_interval(q: F, n=120):
    q = F(q)
    assert q>=1
    z = (q-1)/(q+1)
    if not z: return (F(0),F(0))
    term = z
    acc = F(0)
    for j in range(n):
        acc += 2*term/F(2*j+1)
        term *= z*z
    return (acc, acc + 2*term/(F(2*n+1)*(1-z*z)))

# Analytic inputs: for nonnegative phi, Gamma gap correction >0,
# spacing-three extra prime correlations >=0, and
# 0<=R_phi(log(9/8))<=1/2 since ell<log(11/9)<2log(9/8).
log3_lo = log_interval(F(3))[0]
assert log3_lo > F(109,100)
assert F(5,3)**2<3 and F(7,4)**2>3
# log2 < 7/10, via e^(7/10)>1+x+x²/2+x³/6>2;
# sqrt8 >14/5, hence kappa8=log2/sqrt8 <1/4.
x = F(7,10)
assert 1+x+x*x/2+x*x*x/6 > 2
assert (F(14,5))**2 < 8
coarse_gap = F(5,3)/73*(10*F(109,100)-32*F(1,4))
assert coarse_gap==F(29,438) and coarse_gap>F(1,16)
assert (F(11)-4*F(7,4))/22*F(1,16)==F(1,88)

# Schur form identity for two positive selected directions and one
# complementary direction (real slice; complex uses modulus squares).
a, b = sp.symbols('a b', positive=True)
x0,x1,t0,t1,z=sp.symbols('x0 x1 t0 t1 z',real=True)
D=sp.Symbol('D', real=True)
form = a*x0**2+b*x1**2+2*t0*x0*z+2*t1*x1*z+D*z**2
completed = a*(x0+t0*z/a)**2+b*(x1+t1*z/b)**2+(D-t0**2/a-t1**2/b)*z**2
assert sp.simplify(form-completed)==0
assert sp.det(sp.Matrix([[1,2],[2,1]]))==-3

print('PASS: full spacing-window census through next spacing-two prime 11')
print('Threshold sequence:',', '.join(name for name,_ in thresholds))
print('PASS: two-Mellin-pole annihilator and original physical metric')
print('PASS: outer prime-31 B3 compression',B3)
print('PASS: split increment coefficient',split_coefficient,'* kappa31 * R_phi(delta31)')
print('PASS: nonnegative-profile uniform splitting >',coarse_gap,'> 1/16 through width < log(11/9)')
print('PASS: source-form Schur identity (real slice) and cross-block countermodel')
print('NO proof of wider-width restricted positivity, infinite complement coercivity, or RH')
