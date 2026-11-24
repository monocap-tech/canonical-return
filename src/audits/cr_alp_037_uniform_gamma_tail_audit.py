#!/usr/bin/env python3
"""CR-ALP-037 exact arithmetic for full infinite-sine original Gamma tails.

Continuum resolvent / pole projector / translation identities are analytical
inputs, NOT Lean theorems. This script uses only exact Fraction comparisons for
all asserted bounds and is intentionally nonpromoting. No 15k Schur certificate.
"""
from fractions import Fraction as F
from math import factorial, isqrt

A = F(11, 10)
N = 15000
M = 10000


def exp_taylor_lower(x, terms):
    return sum((x**k / factorial(k) for k in range(terms)), F(0))

# All exponent comparisons are proved with a FINITE exact Taylor lower bound.
assert exp_taylor_lower(F(11, 2), 9) > 200
assert exp_taylor_lower(A, 12) > 3
assert exp_taylor_lower(4*A, 20) > 80
print('PASS exponential enclosures: e^(11/2)>200, e^a>3, e^(4a)>80')

# Interior source columns: exact m=0 cancellation; the first omitted resolvent
# has beta=5/2 and 1+e^(-2a beta)<201/200. The pole-subtracted coefficient
# satisfies |gamma_jk-t_k c_j|<201*k/(200*j^2) for all j>N.
k2sum = sum(k*k for k in range(3, 13))
assert k2sum == 645
interior_gamma_sq = F(201, 200)**2 * k2sum / (3*N**3)
assert interior_gamma_sq < F(1, 120000)**2
# Carry forward the independent source-indexed original PRIME bound from 036.
# This is a theorem-use check, not independent recomputation of CR-ALP-036.
assert F(1, 120000) + F(1, 18000) == F(23, 360000)
print('PASS gamma interior 10-column L2 operator norm <1/120000')
print('PASS combined original gamma + signed seven-prime interior <23/360000')

# Fraction enclosures of genuine original pole-anchor gamma coefficients.
# 3.14<pi<22/7; 26/25<sqrt(11/10). Taylor hyperbolic enclosures.
pi_lo, pi_hi, sqrt_a_lo = F(157, 50), F(22, 7), F(26, 25)
assert pi_lo < pi_hi and sqrt_a_lo**2 < A
x = A/2

def cosh_upper(x):
    terms = sum((x**(2*j)/factorial(2*j) for j in range(5)), F(0))
    return terms+(x**10/F(factorial(10)))/(1-x*x/F(11*12))

def sinh_upper(x):
    terms = sum((x**(2*j+1)/factorial(2*j+1) for j in range(5)), F(0))
    return terms+(x**11/F(factorial(11)))/(1-x*x/F(12*13))

sinh_a_lower = sum((A**(2*j+1)/factorial(2*j+1) for j in range(4)), F(0))
assert cosh_upper(x)<F(6, 5) and sinh_upper(x)<F(3, 5)
assert 2*A < A+sinh_a_lower
assert sinh_a_lower > A

# Gamma anchor u_k = e_k-gamma_k b_k. Bound gamma's from original exact
# sine-pole overlap 2nu_k*b_k(a)/(sqrt(a)(nu_k^2+1/4)) and pole norm.
for kap, gam_bound in [(1,F(3,5)), (2,F(5,3))]:
    nu_lo = pi_lo*kap/(2*A)
    nu_hi = pi_hi*kap/(2*A)
    b_up = cosh_upper(x) if kap==1 else sinh_upper(x)
    B_lo = A+sinh_a_lower if kap==1 else sinh_a_lower-A
    gamma_interval_upper = 2*nu_hi*b_up/(sqrt_a_lo*(nu_lo**2+F(1,4))*B_lo)
    assert gamma_interval_upper < gam_bound
    print(f'PASS original anchor gamma_{kap} < {gam_bound} (rational enclosing value {float(gamma_interval_upper):.9f})')

# For a pole profile b=cosh(x/2),sinh(x/2), let B0>=||b||_infty,
# B1>=||b'||_infty, b_0=b(-a), I_beta=integral_0^(2a) e^-beta*t b(-a+t)dt.
# |I_beta-b_0/beta| <= B1/beta^2+B0*exp(-2a beta)/beta.
# Summing m>=0 with beta_m=2m+1/2 gives S_E <= (9/2)B1+(160/237)B0.
# At M=10000, subtract t_M*b where t_M=sum_{m<=M} 1/beta_m.
# Exact J_beta(b) sine coefficient identity then bounds the low-beta
# subtraction and the remaining m>M gamma kernels.
exp_harmonic_ceiling=F(160,237)   # 2*(1/3)/(1-1/80)
S_beta = F((M+1)*(2*M+1),2) + F(1,4)*(2+F(M,2))
nu_N_lower = F(3*N,2)/A
assert nu_N_lower > 2*M
assert F(11,7)/A*N/(2*M) == F(15,14)

for kap, B0, B1, gam_bound, target in [
    (1,F(6,5),F(3,10),F(3,5),F(1,20)),
    (2,F(3,5),F(3,5),F(5,3),F(1,9)),
]:
    SE = F(9,2)*B1+exp_harmonic_ceiling*B0
    C0 = 2*SE + 2*B0*S_beta/(nu_N_lower**2) + F(27,8)*B0+3*B0/F(14)
    C1 = 3*B0
    # forall j>N: |r_{j,b}| <=(11/(15j))*(C0+C1*log(j/N));
    # sum_{j>N}(C0+C1*log(j/N))^2/j^2 <=
    #    (C0^2+2C0*C1+2C1^2)/N via decreasing integrand.
    b_tail_sq = F(11,15)**2*(C0*C0+2*C0*C1+2*C1*C1)/N
    # For e_k (the sine anchor) the pole-subtracted gamma tail is
    # bounded by (201/200)*k/sqrt(3*N^3), always <1/10^6.
    sine_tail_sq = F(201,200)**2*kap*kap/(3*N**3)
    assert sine_tail_sq < F(1,1000000)**2
    # Triangle inequality without machine sqrt: exact rational bound.
    assert gam_bound**2 * b_tail_sq < (target-F(1,1000000))**2
    print(f'PASS original u_{kap} gamma Q_15000 residual <{target} (source/interval arithmetic)')
    print('  pole correction constants C0,C1:',float(C0),float(C1))
    print('  squared b-only tail bound:',float(b_tail_sq))

# Independent exact algebra check for endpoint-renormalized Gamma kernel:
# c_j =2*b(-a)*nu/(sqrt(a)*(nu^2+1/4)), I_beta=b(-a)/beta+E.
# We use sqrt(a)=1 in this purely formal rational regression.
for nu in (F(7),F(19),F(15001)):
    for beta in (F(1,2),F(5,2),F(41,2)):
        b0,E=F(-3,5),F(2,17)
        c=2*b0*nu/(nu*nu+F(1,4))
        I=b0/beta+E
        lhs=2*nu*nu/(beta*(beta*beta+nu*nu))*c-2*nu/(beta*beta+nu*nu)*I-c/beta
        rhs=-2*nu/(nu*nu+beta*beta)*(E+b0*(beta*beta+F(1,4))/(beta*(nu*nu+F(1,4))))
        assert lhs==rhs
print('PASS exact endpoint-renormalized J_beta(b)-c_j/beta algebra')

# Complete original seven-prime boundary-anchor tails: bounded variation.
# For zero-extended u_k=e_k-gamma_k*b_k, Var_R(e_k)=2*k/sqrt(a).
# Var_R(b_e)=2(cosh(a/2)-1)+2*cosh(a/2);
# Var_R(b_o)=2*sinh(a/2)+2*sinh(a/2).
# Using sqrt(a)>1, gamma_1<3/5, gamma_2<5/3,
# cosh(a/2)<6/5 and sinh(a/2)<3/5 gives Var_R(u_1)<4,
# Var_R(u_2)<8. Signed prime translations have total BV <=2*S_a*Var(u_k).
assert F(2)+F(3,5)*F(2,5)+2*F(3,5)*F(6,5) < F(4)
assert F(4)+F(5,3)*F(6,5)+2*F(5,3)*F(3,5) == F(8)
# Strict bounds for the underlying gamma/hyperbolics make Var(u_2)<8.
S_prime = F(177,50)
prime_anchor_bounds=[(1,F(4),F(1,5)),(2,F(8),F(2,5))]
for kap, V, target in prime_anchor_bounds:
    # IBP on I with zero-extended BV boundary contributions:
    # |<e_j,Pu>| <= 2*S_a*V/(sqrt(a)*nu_j)
    #                   < (2a/3)*(2*S_a*V)/j.
    squared = ((2*A/F(3))*2*S_prime*V)**2/F(N)
    assert squared < target**2
    print(f'PASS original seven-prime u_{kap} high-tail residual <{target}')
assert F(1,20)+F(1,5)==F(1,4)
assert F(1,9)+F(2,5)==F(23,45)
print('PASS complete signed original u_1 residual <1/4, u_2 residual <23/45')

# Original seven-prime sign and authentic support remain source constraints.
weights={2:2,3:3,4:2,5:5,7:7,8:2,9:3}
assert set(weights)=={2,3,4,5,7,8,9}
print('CR-ALP-037 PASS exact rational infinite-sine Gamma tail bounds:')
print('  interiors (w_3..w_12) <1/120000 (physical operator norm),')
print('  anchors u_1 <1/20, u_2 <1/9 (individual source columns).')
print('OPEN: full physical 12-column combined operator (mass-conditioning),')
print('      12..15000 signed shell, 15000 finite Schur, Lean and RH.')
