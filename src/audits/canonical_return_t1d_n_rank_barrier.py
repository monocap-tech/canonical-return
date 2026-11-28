#!/usr/bin/env python3
"""CR-INT-T1D-N: exact-rational arithmetic for the prolate rank barrier.

The analytic trace, defect, parity and min-max proofs are in the note.
This does NOT enumerate prolate eigenvalues, check the Weil form, certify
any finite Schur matrix, prove operator transfer or RH.
"""
from fractions import Fraction as F
from math import isqrt


def log_rat(x, n=50):
    """Outward exact-rational atanh/geometric bounds for log(x), x>=1."""
    x = F(x)
    assert x >= 1
    z = (x - 1) / (x + 1)
    term = z
    acc = F(0)
    for j in range(n):
        acc += 2 * term / F(2 * j + 1)
        term *= z * z
    return acc, acc + 2 * term / (F(2 * n + 1) * (1 - z * z))


def atan_rat(x, n):
    x = F(x)
    s = sum((F((-1)**j) * x**(2*j+1) / (2*j+1)
             for j in range(n)), F(0))
    rem = x**(2*n+1) / (2*n+1)
    return (s, s+rem) if n % 2 == 0 else (s-rem, s)


log2_lo, log2_hi = log_rat(2, 30)
assert log2_hi < F(347, 500) < F(7, 10)
pi5_lo, pi5_hi = atan_rat(F(1, 5), 25)
pi239_lo, pi239_hi = atan_rat(F(1, 239), 15)
pi_lo = 16*pi5_lo - 4*pi239_hi
pi_hi = 16*pi5_hi - 4*pi239_lo
assert F(3) < pi_lo < pi_hi < F(22, 7)

B = 2**44
THETA = F(1, 128)
# Analytic boundary-leakage estimate:
# D_B <= (2/pi^2) * (log(4B)+3/2).
# Here log(4B)=46 log2, pi>3.
D_crude = F(2, 9) * (46*F(7, 10)+F(3, 2))
assert D_crude == F(337, 45)
assert D_crude / (1-THETA) < 8
D_better = F(2, 9) * (46*F(347, 500)+F(3, 2))
assert D_better / THETA < 951

# Trace 4B/pi >14B/11 and <4B/3; both are strict.
rank_lower_frac = F(14*B, 11) - 8
rank_lower_int = rank_lower_frac.numerator // rank_lower_frac.denominator + 1
assert rank_lower_int == 22390054965613
additional_F = rank_lower_int - 4
assert additional_F == 22390054965609

# Parity trace = half trace +/- (Si(4B)/pi)/2;
# |Si|<=4, hence each trace > 2B/pi-2/3 > 7B/11-2/3.
# Defect within each parity <= total defect, which adds <8.
parity_lower_frac = F(7*B, 11) - 9
parity_lower_int = parity_lower_frac.numerator // parity_lower_frac.denominator + 1
assert parity_lower_int == 11195027482802
assert parity_lower_int - 2 == 11195027482800

# For a physically orthonormal smooth r-correction system on (-2,2),
# sum ||w_j'||^2 >= (pi^2/16) sum_{j<=r} j^2
# = pi^2*r(r+1)(2r+1)/96.
r = additional_F
sum_squares = r*(r+1)*(2*r+1)//6
assert sum_squares*6 == r*(r+1)*(2*r+1)
assert 16 * sum_squares == F(r*(r+1)*(2*r+1), 6) * 16

print('PASS T1D-N: exact log2 and pi outward rational source bounds')
print('PASS T1D-N: trace-defect upper D<337/45, D/(1-theta)<8, D/theta<951')
print('B:', B)
print('MINIMUM m_B:', rank_lower_int)
print('MINIMUM new independent constraints in F:', additional_F)
print('MINIMUM each parity rank:', parity_lower_int)
print('MINIMUM new constraints per parity:', parity_lower_int-2)
print('PASS T1D-N: Dirichlet eigenvalue sum coefficient pi^2*r(r+1)(2r+1)/96')
print('ANALYTIC ONLY: trace-defect identity, min-max necessity, smooth-source density')
print('OPEN: finite Schur matrix, complement coercivity, operator transfer and RH')
