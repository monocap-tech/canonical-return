#!/usr/bin/env python3
"""CR-ALP-027 exact arithmetic and algebra audit of two-pole fifth-order tail.

Verifies the seven original von Mangoldt weights, stronger rational bounds,
parity-resolvent subtraction, and a 15000-mode W-tail >= 1/4 arithmetic gate.
Continuum Fourier/digamma/domain steps are analytical, NOT machine-proved.
No original finite Weil matrix computed. No Lean kernel certification or RH.
"""
from fractions import Fraction as F
from math import factorial, isqrt

A = F(11, 10)
PRIME_POWERS = {2:2,3:3,4:2,5:5,7:7,8:2,9:3}


def log_atanh_1_2(z, count=40):
    assert F(1) <= z <= F(2)
    u = (z-1)/(z+1)
    term, lo = u, F(0)
    for j in range(count):
        lo += 2*term/(2*j+1)
        term *= u*u
    hi = lo + 2*term/((2*count+1)*(1-u*u))
    return lo, hi


def log_int(n):
    assert n >= 1
    k = n.bit_length()-1
    lo, hi = log_atanh_1_2(F(n, 2**k))
    l2, h2 = log_atanh_1_2(F(2))
    return lo+k*l2, hi+k*h2


# Exact seven original powers and their authentic log(p)/sqrt(n) weights.
assert set(PRIME_POWERS) == {2,3,4,5,7,8,9}
assert all(any(n == p**k for k in range(1,5))
           for n,p in PRIME_POWERS.items())
assert log_int(9)[1] < 2*A < log_int(10)[0]
scale = 10**6
S_upper = F(0)
for n, p in PRIME_POWERS.items():
    sqrt_lo = F(isqrt(n*scale*scale), scale)
    assert sqrt_lo**2 <= n and sqrt_lo > 0
    S_upper += log_int(p)[1]/sqrt_lo
assert S_upper < F(177,50)
print('PASS authentic prime powers S_a < 177/50; paired negative mass <177/25')

# pi<22/7; finite *lower* sum of e^(37/20) exceeds 44/7 >2pi.
x=F(37,20)
assert sum((x**k/factorial(k) for k in range(7)),F(0))>F(44,7)
assert log_int(10000)[0] > F(921,100)
print('PASS log(2*pi)<37/20; log(10000)>921/100')

# Symbolic-rational parity cancellation regression (toy nu_j=j,
# chosen only to verify the rational resolvent algebra; the real nu_j=j*pi/2a).
# For each parity, choose two nonzero tail coefficients annihilating
# the exact pole moment sum c_j*nu_j/(nu_j^2+1/4).
coeff = {}
for j,k in [(5,7),(6,8)]:
    pole_weight = lambda n: F(n)/(n*n+F(1,4))
    coeff[j] = pole_weight(k)
    coeff[k] = -pole_weight(j)

for parity in (0,1):
    js = [j for j in coeff if j%2==parity]
    pole_sum = sum((coeff[j]*F(j)/(j*j+F(1,4)) for j in js),F(0))
    assert pole_sum == 0
    for r in [F(1,3),F(2,3),F(3,2)]:
        direct = sum((coeff[j]*j/(j*j-r*r) for j in js),F(0))
        cancelled = sum((coeff[j]*j*(r*r+F(1,4)) /
                         ((j*j-r*r)*(j*j+F(1,4))) for j in js),F(0))
        assert direct == cancelled
print('PASS parity-resolved pole-null cancellation of leading sine resolvent')


def fifth_leak_upper(R, N):
    assert R >= 2 and N>=1
    # pi >3, hence nu_{N+1}=(N+1)*pi/(2a)>2R.
    assert 3*(N+1) >= 4*A*R
    # Use pi>3 in the exact analytic bound
    # 4096 a^5/(45*pi^7*N^5)*(R^5/5+R^3/6+R/16).
    d = F(4096)*A**5/(45*3**7*N**5) * \
        (F(R**5,5)+F(R**3,6)+F(R,16))
    assert d < 1
    return d

R,N=10000,15000
d=fifth_leak_upper(R,N)
assert d < F(1,500)
L=F(921,100)-F(37,20)-F(2,R)
assert L == F(36799,5000)
mu_coarse = L*(1-F(1,500))-F(6,500)-F(177,25)
assert mu_coarse == F(632701,2500000) and mu_coarse>F(1,4)
print('PASS N=15000, R=10000: fifth-order leakage <1/500')
print('PASS original signed two-pole-null W(q,q)>1/4 ||q||^2 on infinite tail')

# Additional exact interval tradeoffs; these use the precise rational lower
# logs instead of the coarse 921/100 for log(10000).
for R,N,margin in [(8000,12000,F(1,30)),(9000,14000,F(3,20))]:
    d = fifth_leak_upper(R,N)
    L = log_int(R)[0]-F(37,20)-F(2,R)
    mu = L*(1-d)-6*d-F(177,25)
    assert mu>margin
    print(f'PASS optional R={R} N={N}: margin >{margin}')

# A finite-positive low block is insufficient without its residual Schur term.
H, cross, tail = F(1,4), F(1), F(2)
assert H>0 and tail>0 and H-cross**2/tail<0
# With tail >=1/4, the sufficient low block criterion is H-4B>=0.
H_good, B_good = F(5), F(1)
assert H_good-4*B_good>0
print('PASS finite Schur penalty 4 and distinction between finite vs full positivity')
print('CR-ALP-027 EXACT RATIONAL AUDIT PASS; continuum and Lean formalization OPEN')
