#!/usr/bin/env python3
"""CR-ALP-007: exact rational bounds for three translated original Weil ramps.

External analytic premises: standard Euler-constant enclosure by H_N-log N,
Archimedes' 223/71 < pi < 22/7, and elementary antiderivatives.
ALL certification inequalities use Fraction, never floating point.
No Lean compilation, theta-lift triple, operator completion or RH theorem.
"""
from fractions import Fraction as Q


def primes_through(n):
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, int(n ** .5) + 1):
        if sieve[p]:
            for j in range(p * p, n + 1, p):
                sieve[j] = False
    return [i for i in range(2, n + 1) if sieve[i]]


def log_atanh(q, N=18):
    """Rigorous rational bounds for log q, with q in [1,2]."""
    assert Q(1) <= q <= Q(2)
    z = (q - 1) / (q + 1)
    acc, term = Q(0), z
    for j in range(N):
        acc += 2 * term / (2 * j + 1)
        term *= z * z
    err = 2 * term / ((2 * N + 1) * (1 - z * z))
    return acc, acc + err


def log_int(n, N=18):
    assert n >= 1
    m = n.bit_length() - 1
    a, b = log_atanh(Q(n, 2 ** m), N)
    c, d = log_atanh(Q(2), N)
    return a + m * c, b + m * d


def log_ratio(a, b=1, N=18):
    l, u = log_int(a, N)
    l2, u2 = log_int(b, N)
    return l - u2, u - l2


def enclose(terms, N=18):
    lo = hi = Q(0)
    for coeff, a, b in terms:
        ell, upper = log_ratio(a, b, N)
        if coeff >= 0:
            lo += coeff * ell
            hi += coeff * upper
        else:
            lo += coeff * upper
            hi += coeff * ell
    return lo, hi


def coefficients(limit, weight):
    out = {}
    for p in primes_through(limit):
        n, cp = p, Q(0)
        while n <= limit:
            cp += weight(n)
            n *= p
        if cp:
            out[p] = cp
    return out


def first_weight(n):
    if n <= 9:
        return Q(3 * (n - 1), n)
    return Q(81 - n, 3 * n)


def second_weight(n):
    if n <= 9:
        return Q(0)
    if n <= 81:
        return Q(1) - Q(9, n)
    return Q(81, n) - Q(1, 9)


# H01=H12, source separation log 9; complete two-pole/prime/Gamma formula.
c01 = coefficients(81, first_weight)
assert len(c01) == 22
assert c01[2] == Q(533, 64) and c01[3] == Q(16, 3)
terms01 = [(Q(27), 9, 1)]
terms01 += [(-c, p, 1) for p, c in c01.items()]
terms01 += [(Q(-3), 5, 1), (Q(-27, 2), 50, 41), (Q(1, 6), 82, 81)]
lo01, hi01 = enclose(terms01)
assert Q(-132426, 1000000) < lo01 < hi01 < Q(-132425, 1000000)

# H02: separation 2 log 9; prime powers in 9<n<=729.
c02 = coefficients(729, second_weight)
assert len(c02) == 129
assert c02[2] == Q(4285, 1536) and c02[3] == Q(16, 9)
terms02 = [(Q(72), 9, 1)]
terms02 += [(-c, p, 1) for p, c in c02.items()]
terms02 += [
    (Q(-1, 2), 82, 81),
    (Q(1, 18), 6643, 6642),
    (Q(9, 2), 50, 41),
    (Q(81, 2), 3650, 3731),
]
lo02, hi02 = enclose(terms02)
assert Q(-1978, 1000000) < lo02 < hi02 < Q(-1977, 1000000)

# H00=H11=H22: complete original signed Weil diagonal for cutoff X=9.
cdiag = coefficients(9, lambda n: 2 * (Q(9, n) - 1))
assert cdiag == {2: Q(39,4), 3: Q(4), 5: Q(8,5), 7: Q(4,7)}
termsd = [(Q(16), 9, 1)]
termsd += [(-c, p, 1) for p, c in cdiag.items()]
termsd += [(Q(10), 5, 1), (Q(-16), 4, 1)]
lod, hid = enclose(termsd)

# H_100 - log 100 - (2*100)^-1 < Euler gamma
# < H_100 - log 100 - (2*101)^-1.
N = 100
h = sum((Q(1,i) for i in range(1,N+1)), Q(0))
lN, uN = log_int(N)
glo, ghi = h-uN-Q(1,2*N), h-lN-Q(1,2*(N+1))
lod -= 8*ghi
hid -= 8*glo

# Archimedes' rational pi enclosure.
pilo, _ = log_ratio(223,71)
_, pihi = log_ratio(22,7)
lod -= 8*pihi
hid -= 8*pilo
assert Q(450,1000) < lod < hid < Q(456,1000)

# Physical 3x3 metric M=8 I. The signed H is symmetric Toeplitz with
# diagonal d and off-diagonal c1, c2 (at distances L and 2L).
assert hi01 < 0 and hi02 < 0
margin = Q(450,1000) - 2*Q(132426,1000000)
assert margin > Q(185,1000)
assert margin / 8 > Q(23,1000)

# Cyclic Bargmann-type invariant Xi=c1^2*c2/d^3 <0.
chi_abs_lo = Q(132425,1000000)**2 * Q(1977,1000000) / Q(456,1000)**3
chi_abs_hi = Q(132426,1000000)**2 * Q(1978,1000000) / Q(450,1000)**3
assert Q(36,100000) < chi_abs_lo < chi_abs_hi < Q(39,100000)

for name,low,high in [('H01=H12',lo01,hi01),
                      ('H02',lo02,hi02),('diagonal',lod,hid)]:
    print(f'{name}: ({float(low):.12f},{float(high):.12f}) [display only]')

print('CR-ALP-007 PURE-RATIONAL AUDIT PASS')
print('Three physically orthogonal ramp sources, M=8 I_3; signed H>0;')
print('three strictly negative mixed entries and a negative cyclic holonomy.')
print('No theta-lift operator-instance, Lean, semilocal or RH claim.')
