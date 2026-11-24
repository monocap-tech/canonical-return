#!/usr/bin/env python3
"""CR-ALP-036: original seven-prime cross-tail via exact pole cancellation.

All final bounds use Fraction, integer square roots, rational pi/log/prime
weight envelopes inherited from CR-ALP-027. The analytical identities for
the original shift operator, the infinite sine basis, and Gamma decomposition
are separately stated in the note; symbolic finite regressions test algebra.
No positive signed Weil Schur result, Lean proof, or RH is claimed.
"""
from fractions import Fraction as F
from math import isqrt

A = F(11,10)
S = F(177,50)       # exact audited S_a upper for original seven prime powers
N = 15000
assert N >= 2*12


def inv_sqrt_upper(n, scale=10**12):
    """Strict rational upper bound on 1/sqrt(n), for positive integer n."""
    assert isinstance(n,int) and n > 0
    q = isqrt(n*scale*scale) + 1
    r = F(q,n*scale)
    assert r*r*n > 1
    return r


# The genuine pole-null coefficients are lambda_k =
# x(y^2+1/4)/(y(x^2+1/4)), x=nu_k, y=nu_{kappa}, kappa=1/2.
# Using pi>3 yields lambda_k <= kappa/k + a^2/(9*k*kappa).
# Using pi<22/7 gives nu_k <= (11/7a) k.
# With j>N>=2k, |nu_j^2 - nu_k^2| >= 3nu_j^2/4.
# Prime signed shift exact offdiag U(j,k,t) is (same parity):
# 2( nu_j sin(nu_k t) - nu_k sin(nu_j t) )/[a(nu_j^2-nu_k^2)]
# For w_k=e_k-lambda_k e_kappa, subtract a source-specific multiple
# of the corresponding pole b so the O(1/j) boundary term cancels.
# Then |high_sine_prime_coeff_minus_pole| <= Ak/j^2+Bk/j^3.

bounds=[]
for k in range(3,13):
    kap=1 if k%2 else 2
    lambda_up=F(kap,k)+A*A/F(9*k*kap)
    assert 0 < lambda_up < 1
    x_up=F(11*k,7)/A
    y_up=F(11*kap,7)/A
    Ak=F(8)*S/(3*A)*(x_up+lambda_up*y_up)*(2*A/3)**2
    Bk=((F(8)*S/(3*A)*(x_up**2+lambda_up*y_up**2)
          +S*(1+lambda_up)/(2*A))*(2*A/3)**3)
    # Minkowski + integral test for sum_{j>N} j^{-4},j^{-6}.
    bound=Ak*inv_sqrt_upper(3*N**3)+Bk*inv_sqrt_upper(5*N**5)
    assert bound < F(1,30000)
    bounds.append(bound)
    print(f'k={k:2} lambda<{float(lambda_up):.9f}; '
          f'||Q_N P w_k||<{float(bound):.12g}')

assert bounds[0] < F(1,100000)
assert sum((b*b for b in bounds),F(0)) < F(1,18000)**2
print('PASS: each authentic signed seven-prime interior column <1/30000')
print('PASS: for p=sum_{k=3}^{12} x_k w_k, '
      '||Q_15000 P p|| < ||p||/18000 (physical L2)')

# Algebra-only checks of the source-specific Gamma resolvent cancellation.
# They use rational toy nu's; universal identity follows by clearing the
# rational denominators in the companion mathematical note.
for k in range(3,13):
    kap=1 if k%2 else 2
    x,y=F(13*k,7),F(13*kap,7)
    lam=x*(y*y+F(1,4))/(y*(x*x+F(1,4)))
    for beta in [F(1,2),F(5,2),F(9,2),F(67,2)]:
        direct=x/(beta*beta+x*x)-lam*y/(beta*beta+y*y)
        factor=x*(x*x-y*y)*(beta*beta-F(1,4))/(
            (x*x+F(1,4))*(beta*beta+x*x)*(beta*beta+y*y))
        assert direct==factor
        if beta==F(1,2): assert direct==0
print('PASS: first original Gamma resolvent J_{1/2} has ZERO high-sine '
      'coupling with every genuine interior pole-null w_k')

# A rigorous finite-index resolvent remainder modulus for fixed j,k:
# for M>=2, from |1-sigma e^{-2a beta}|<=2, beta_m>=2m,
#  |sum_{m>=M} J_beta_m(e_j,w_k)|
#    <= nu_j*(nu_k+lambda*nu_kappa)/(12a*(M-1)^3).
# This is explicit for every fixed high index j; not uniform as j->infty.
for k in (3,7,11,12):
    kap=1 if k%2 else 2
    lam=F(kap,k)+A*A/F(9*k*kap)
    j=15001 if k%2 else 15002
    nu_j=F(11*j,7)/A
    nu_k=F(11*k,7)/A
    nu_kap=F(11*kap,7)/A
    M=30000
    rem=nu_j*(nu_k+lam*nu_kap)/(12*A*(M-1)**3)
    assert rem < F(1,100_000_000)
print('PASS: high-index single-column Gamma M=30000 exact tail budget '
      '<1e-8 for four sample source-indexed (j,k) pairs')
# Source-attached boundary anchors u_1,u_2 are NOT in the interior H^1_0
# class. Nevertheless their physical derivative and endpoint budgets admit
# explicit rational bounds: D(u_k)<25 and T(u_k)<4, using
# ||b_e||^2=a+sinh(a), ||b_o||^2=sinh(a)-a,
# <b'_e,b'_e>=||b_o||^2/4 and vice versa, and
# ||gamma_k b|| <= 1 for the exact L2 pole projection.
assert A+2 < 16*A**3/F(6)   # G_e/G_o <16, using sinh a<2.
assert (F(11,7)/A+F(1,2))**2 < 4
assert (F(22,7)/A+2)**2 < 25
assert 2*F(6,5)**2/(2*A) < 4
assert 2*F(3,5)**2/(A**3/F(6)) < 4
M=30000
# For J_beta(v,v) <= 4 D(v)/beta^3 +4 T(v)/beta^2,
# and J_beta(e_j,e_j) <=2nu_j^2/beta^3, Cauchy gives
# |J_beta(e_j,u_k)| <= 2sqrt(2)*nu_j*(sqrt(D)/beta^3
#                                      +sqrt(T)/beta^(5/2)).
# With D<25,T<4, 2sqrt(2)<3, and
# sum_{m>=M} beta_m^-3 <= (M^-3+1/(2M^2))/8;
# sum_{m>=M} beta_m^-5/2 <= (M^-2+M^-1)/800 when M=30000,
# since 1/sqrt(2M)<1/200. This is a PER-INDEX modulus, not a
# global infinite-sine residual estimate.
assert 2*M > 200**2
beta_3=F(1,M**3)/8+F(1,16*M**2)
beta_5_2=(F(1,M**2)+F(1,M))/800
for j in (15001,15002):
    nu_j_upper=F(11*j,7)/A
    anchor_gamma_remaining=3*nu_j_upper*(5*beta_3+2*beta_5_2)
    assert anchor_gamma_remaining < F(1,100)
print('PASS: source boundary-anchor derivative D<25, trace T<4, and '
      'rigorous Gamma high-index tail modulus <1/100 at j=15001/15002')
print('CR-ALP-036 EXACT RATIONAL/ALGEBRAIC AUDIT PASS; full Schur OPEN')
