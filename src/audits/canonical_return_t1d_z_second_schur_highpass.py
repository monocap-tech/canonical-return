#!/usr/bin/env python3
"""CR-INT-T1D-Z: rational second-stage Schur-norm/high-pass audit.

Strictly finite exact-fraction checks. Imports T1D-E/F/Q/U/W/X/Y
source-level analytic estimates and T1D-M's Gamma/prolate leakage lemma.
No full-complement positivity, no independent analytic/formal verification, no RH.
"""
from fractions import Fraction as F

# Rebuild the inherited Y determinant lower certificates from the
# source-level W/X/Y bounds (all are strict, not rounded data).
root2_upper = F(283,200)
assert root2_upper**2 > 2
D = {}
for parity, qyy, qzz, cross, ryy, rzz, ryz, certificate in (
    ('even',F(89,100),F(1159,1000),root2_upper*F(69,100),F(4,1575),F(729,70000),F(9,1750),F(1,20)),
    ('odd',F(118,100),F(1069,1000),root2_upper*F(72,100),F(9,4375),F(529,70000),F(69,17500),F(1,5)),
):
    diag_y,diag_z=qyy-ryy,qzz-rzz
    det=diag_y*diag_z-(cross+ryz)**2
    assert diag_y>0 and diag_z>0 and det>certificate
    D[parity]=det
    print(f'PASS Y imported {parity} Schur determinant lower: {float(det):.9f} > {certificate}')

# W (complete Gamma central two-space upper), X (satellite wing prime
# plus off-diagonal Gamma) yield full q, hence full s, trace bounds.
trace_even=2*F(31,20)+F(29,100)+F(1,1000) # 3.391
trace_odd=2*F(9,5)-F(9,100)+F(1,1000)    # 3.511
assert trace_even < F(7,2) and trace_odd<F(18,5)
# Since each 2x2 parity S positive, lambda_min >= det/trace.
assert F(1,20)/F(7,2)==F(1,70)
assert F(1,5)/F(18,5)==F(1,18)
print('CERTIFIED analytic inputs: lambda_min(S_even)>1/70, lambda_min(S_odd)>1/18')

# T1D-F independent (no T1D-I) one-sided beta_sigma<112, mu>7.
# T1D-F pointwise Gamma multiplier bound implies
#   ||Gamma f|| <= 40||f||+||f_prime||.
# T1D-U central source ||A P_J||<7/2, T1D-Q global ||A||<203/25.
# Central/satellite normalized H1 ratios: <38 even, <46 odd (W/X).
# T1D-W selected-to-central cross: <2/15 even, <3/25 odd.
# T1D-Y selected-to-satellite cross: <27/100 even, <23/100 odd.
for parity,r,cross_c,cross_s,cap_c,cap_s in (
    ('even',38,F(2,15),F(27,100),84,91),
    ('odd',46,F(3,25),F(23,100),92,98),
):
    exact_c=40+r+F(7,2)+cross_c*F(112,7)
    exact_s=40+r+F(203,25)+cross_s*F(112,7)
    assert exact_c<cap_c and exact_s<cap_s,(parity,exact_c,exact_s)
    print(f'PASS {parity} source representers: ||g_center||<{cap_c}, '
          f'||g_satellite||<{cap_s} (rational budgets {exact_c}, {exact_s})')

rank_even=70*(84**2+91**2)
rank_odd=18*(92**2+98**2)
assert rank_even == 1_073_590 and rank_odd == 325_224
rank_bound=max(rank_even,rank_odd)
print('CERTIFIED new rank<=4 penalty operator bound:',rank_bound,
      '(even:',rank_even,'odd:',rank_odd,')')

# All source operators remain REAL reflection-invariant, so the even/odd
# penalty norms combine by max, NOT by sum. Extra four constraints
# b_i(w)=0 annihilate the rank<=4 penalty exactly.
# Check the rational lower bound log2>693/1000 via atanh series.
z=F(1,3)
log2_lower=sum((2*z**(2*j+1)/F(2*j+1) for j in range(12)),F(0))
assert log2_lower>F(693,1000)

def low_weil(k):
    # T1D-M/Q high Gamma bound (for B=2^k, low band theta=1/128),
    # signed *joint* 24-prime operator norm bound 203/25.
    return F(127,128)*(k*F(693,1000)-F(49,10))-F(6,128)-F(203,25)

coarse_R=F(1785) # independent source-level T1D-F, NOT T1D-I
with_extra_conditions=low_weil(2700)-coarse_R
without_extra_conditions=low_weil(1_570_000)-coarse_R-rank_bound
assert with_extra_conditions == F(187051,3200)
assert without_extra_conditions == F(26380117,6400)
assert with_extra_conditions>58 and without_extra_conditions>4000
print('PASS k=2700 with four exact b_i-null constraints: second Schur >',
      with_extra_conditions,'>',58)
print('PASS k=1570000 without b_i-null constraints: second Schur >',
      without_extra_conditions,'>',4000)
print('ANALYTICAL INPUTS: full signed Weil Gamma/prolate lemma; F one-sided beta<112 and R<1785; W/X/Y central-satellite Grams; E source representer; Q/U prime norms')
print('OPEN: low-mode second Schur matrix, uniform positivity on all W, top compact spectrum, RH')
