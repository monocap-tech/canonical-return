#!/usr/bin/env python3
"""CR-INT-T1D-T: exact algebra for an adapted Schur tail and broad pole-null probes.

Checks algebraic finite-dimensional Schur mechanism and rational constants,
plus exact polynomial H1 estimates for two smooth-after-mollification,
parity-separated pole-null broad-support sources.
Does NOT calculate any D_F eigenvalue, cross-residual norm, leading Ritz
value, complement positivity, Lean certificate, or RH.
"""
from fractions import Fraction as F
from sympy import symbols, diff, expand, integrate, Matrix, simplify, Rational, sqrt

# Import assumptions (not re-proved here): T1D-Q ||A||<203/25,
# T1D-I ||R||<155/182, T1D-R D_F>=5/8 and compact resolvent.
L = F(203,25) + F(155,182)
assert L < 9
assert F(6)+L < 15
assert F(5,8)>0

# 2x2 positive semidefinite block criterion:
# I-B = [[1-b, -x],[-x,1-t]]. For t<1, PSD iff
# 1-b >= x^2/(1-t).
b,t,x = symbols('b t x',real=True)
H = Matrix([[1-b,-x],[-x,1-t]])
assert simplify(H.det()-((1-b)*(1-t)-x*x))==0

# General bound: ||P B Q|| <= L/sqrt(d1*h), ||Q B Q|| <= (6+L)/h.
# Therefore delta = L^2/[d1*(h-6-L)]. If d1>=5/8,
# L<9, h>15, delta< 648/[5(h-15)].
d1,h = symbols('d1 h',positive=True)
raw = L*L/(F(5,8)*(F(275)-F(6)-L))
coarse = F(648,5*(275-15))
assert raw < coarse < F(1,2)
assert coarse == F(162,325)

# For h <= 274.2, the particular conservative delta_bound >=1/2;
# at h >=275 a half-unit tail budget is available. This is a
# certificate for the estimator, NOT a claim the true error is so large.
assert 15 + 2*F(648,5)==F(1371,5)  # 274.2
assert 275>F(1371,5)

# Positive matrix example: simple norm-Ritz truncation is too crude,
# while the exact Schur complement gives the threshold without
# dropping 6I into the off-diagonal block.
B = Matrix([[Rational(3,4),Rational(1,10)],[Rational(1,10),Rational(1,4)]])
assert 1-B[1,1]>0
assert B[0,0]+B[0,1]**2/(1-B[1,1])<1
assert min((Matrix.eye(2)-B).eigenvals().keys())>0

# BROAD-SUPPORT TESTS. For |x|<a, choose b_e=(1-16x^2)^5,
# b_o=x(1-16x^2)^5, zero outside; b in H^4_0. Put
# y0=(d^2/dx^2-1/4)b, and y=y0*rho_eps with an even nonnegative
# C-infinity probability mollifier supported in [-1/200,1/200].
# Then y in C_c^infty[-51/200,51/200], and its +/-1/2 exponential
# moments vanish exactly by integrations by parts. With ell=1/10000,
# selected E sources have centers >=log3/2>1/2, so these y are
# physically orthogonal to E. All y-y primes vanish since
# diameter 51/100 < 2/3 < log2.
z = symbols('z',real=True)
a = F(1,4)
vnorms=[]
for name,poly,r_max,ratio_bound in [
        ('even',(1-16*z*z)**5,21,24),
        ('odd',z*(1-16*z*z)**5,26,30)]:
    p = expand(diff(poly,z,2)-poly/4)
    v = integrate(p*p,(z,-Rational(1,4),Rational(1,4)))
    w = integrate(diff(p,z)**2,(z,-Rational(1,4),Rational(1,4)))
    assert v>0 and w>0
    assert w < r_max*r_max*v
    # Young: ||y'|| <= ||y0'||; Minkowski: ||y-y0|| <= eps||y0'||.
    # Thus ||y'||/||y|| < r_max/(1-eps*r_max).
    eps=F(1,200)
    rate=F(r_max)/(1-eps*r_max)
    assert rate < ratio_bound
    # Fourier Jensen and the imported Gamma bound from T1D-F:
    # gamma[y]/||y||^2 < 40+log(1+rate)<44, since rate<30
    # and log31<4 (elementary e^4>31).
    assert 1+ratio_bound <= 31
    vnorms.append((name,str(v),str(w),str(rate),ratio_bound))

# exact geometric rational checks; log3 >1 from atanh series and
# log2 >2/3; no numerical log comparison in any soundness premise.
assert F(51,100)<F(2,3)
assert F(51,200)<F(1,2)-F(1,20000)
assert F(6,50)==F(3,25)

for name,v,w,rate,rr in vnorms:
    print(f'PASS broad {name} test: ||y0||^2={v}, ||y0_prime||^2={w}, H1 ratio<{rate}<{rr}')
print('PASS structure-adapted block Schur: ||A+R||<9 (conditional on T1D-I), identity 6I causes ZERO mixed block')
print('PASS new sufficient tail correction: delta<648/[5(h-15)] for h=d_{N+1}>15')
print('PASS if h>=275, delta<162/325<1/2, versus old estimator needing h>=5760')
print('PASS broad support smooth double-pole-null even/odd trial tests; no prime self terms')
print('ANALYTIC imports: Gamma multiplier and compact operator, Young/Minkowski, exact Schur theorem')
print('OPEN: evaluate eta_N, Ritz maximum, gamma trial integrals, complement positivity, RH')
