#!/usr/bin/env python3
"""CR-ALP-026 exact Fraction audit: stronger original digamma/prime tail.

Checks arithmetic only. The Euler sum-vs-integral digamma inequality,
Fourier sine leakage theorem and unbounded-operator form statements are
separate human analytic proofs. Not Lean, F41 certification or RH.
"""
from fractions import Fraction as F
from math import isqrt, factorial

A = F(11, 10)
PP = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}


def atanh_log_unit(z, count=40):
    assert F(1) <= z <= F(2)
    u = (z - 1) / (z + 1)
    power, lower = u, F(0)
    for k in range(count):
        lower += 2 * power / (2 * k + 1)
        power *= u * u
    upper = lower + 2 * power / ((2 * count + 1) * (1 - u * u))
    return lower, upper


def log_int(n):
    assert n >= 1
    m = n.bit_length() - 1
    lo, hi = atanh_log_unit(F(n, 2**m))
    l2, h2 = atanh_log_unit(F(2))
    return lo + m * l2, hi + m * h2


# Exact original prime powers, with their genuine Lambda(p^k)=log(p).
assert set(PP) == {2, 3, 4, 5, 7, 8, 9}
assert all(any(n == p**k for k in range(1, 5)) for n, p in PP.items())
assert log_int(9)[1] < 2 * A < log_int(10)[0]
scale = 10**6
S_upper = F(0)
for n, p in PP.items():
    sqrt_low = F(isqrt(n * scale**2), scale)
    assert 0 < sqrt_low and sqrt_low**2 <= n
    S_upper += log_int(p)[1] / sqrt_low
assert S_upper < F(71, 20)
print('PASS: genuine seven prime-power weights; S < 71/20')

# High-frequency Gamma floor: g(r)>= log |r| - log(2*pi) - 2/|r|.
# Source-side continuum input: Re psi(1/4+it)>=log sqrt(1/16+t^2)-1/t.
# Exact elementary comparison log(2*pi)<2: pi<22/7 and e^2>7>44/7.
assert F(44, 7) < 7
exp2_partial = sum((F(2)**k / factorial(k) for k in range(6)), F(0))
assert exp2_partial > 7
# Native lower Gamma floor at all frequencies is >-6, as CR-ALP-024.
assert log_int(3)[0] > 1
assert log_int(10)[0] > F(16, 7)

R = 30_000
L = F(57, 7) - F(2, R)
assert L == F(854993, 105000)
assert F(1) + 4 * F(16, 7) - 2 - F(2, R) == L


def tail_data(N):
    # nu_{N+1}=(N+1)*pi/(2a)>2R, using pi>3.
    assert F(3) * (N + 1) > 4 * A * R
    # delta <= (256*a*R)/(9*pi^3*N) < (256*a*R)/(243*N).
    d = F(256) * A * R / (243 * N)
    assert d < 1
    mu = L * (1 - d) - 6 * d - F(71, 10)
    return d, mu


d500, mu500 = tail_data(500_000)
assert d500 == F(704, 10125)
assert mu500 == F(63181553, 1063125000)
assert mu500 > F(1, 20)
print('PASS: R=30000, N=500000; leakage <', d500,
      '; original signed tail margin >', mu500, '>1/20')

d1000, mu1000 = tail_data(1_000_000)
assert d1000 == F(352, 10125)
assert mu1000 == F(585899089, 1063125000)
assert mu1000 > F(1, 2)
print('PASS: R=30000, N=1000000; original signed tail margin >',
      mu1000, '>1/2')

# The algebraic threshold for positive lower margin at this fixed R.
threshold = (F(256) * A * R / 243) * (L + 6) / (L - F(71, 10))
assert threshold == F(4181740288000, 8868933)
assert threshold < 500_000

# Abstract exact 2x2 Schur witness: finite H>=0 alone cannot certify tail.
H, cross, tail = F(1, 4), F(1), F(2)
assert H > 0 and tail > 0 and H - cross**2 / tail < 0
# New sufficient penalty for N=500k is mu^{-1}=20, not the exact tail inverse.
H_good, residual = F(21), F(1)
assert H_good - 20 * residual > 0
print('PASS: finite Schur logic, size/margin tradeoff, exact rational threshold')
print('CR-ALP-026 EXACT ARITHMETIC PASS; continuum and operator claims unformalized')
