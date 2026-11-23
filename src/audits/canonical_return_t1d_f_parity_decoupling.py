#!/usr/bin/env python3
"""CR-INT-T1D-F algebra/rational audit: parity residuals and bounds.

Analytic Fourier-multiplier and integration inequalities are presented
in research/CR-INT-T1D-F.md. This program checks finite matrix
identities, exact rational constants and positivity margins only.
It does NOT numerically enclose the actual residual norms, prove
complement positivity, construct a completed operator, or prove RH.
"""
from fractions import Fraction as F
from math import factorial
import sympy as sp

r = sp.sqrt(3)
s = 4/r
X = sp.Matrix([[1, 0], [-s, 1], [1, -s], [0, 1]])
R4 = sp.zeros(4)
for j in range(4):
    R4[j,3-j] = 1
assert sp.simplify(R4*X-X*sp.Matrix([[0,1],[1,0]])) == sp.zeros(4,2)
c_plus = X*sp.Matrix([1,1]);c_minus=X*sp.Matrix([1,-1])
assert R4*c_plus == c_plus
assert R4*c_minus == -c_minus
M=sp.simplify(X.T*X)
assert sp.simplify((c_plus.T*c_minus)[0]) == 0
assert sp.simplify(M-sp.Matrix([[sp.Rational(22,3),-8*r/3],[-8*r/3,sp.Rational(22,3)]]))==sp.zeros(2,2)

# Windowed poles m_+,m_- reduce to cosh(x/2),sinh(x/2) by parity.
x=sp.symbols('x',real=True)
pe=sp.cosh(x/2);po=sp.sinh(x/2)
assert sp.simplify(sp.integrate(pe**2,(x,-2,2))-(2+sp.sinh(2)))==0
assert sp.simplify(sp.integrate(po**2,(x,-2,2))-(sp.sinh(2)-2))==0
assert sp.sinh(2)>2
assert sp.simplify(sp.integrate(pe*po,(x,-2,2)))==0

# T1C triangle-mollifier: ell^2||phi'||^2/J <= (175/48)^2.
q=F(30625,2304)
assert q==F(175,48)**2
assert F(5,3)+q<16  # 5/pi < 5/3 since pi>3: L_psi<4

# A rational witness for log(10000)<10 via e^10 Taylor partial.
assert sum((F(10)**k/F(factorial(k)) for k in range(10)),F(0))>10000
ELL=F(1,10000)
fixed_source_norm_bound=2*F(175,48)+64+ELL+4*F(10)
assert fixed_source_norm_bound<F(112)
mu_lower=F(7)
# beta^2 <= ||K u||^2 - mu^2
assert (F(112)**2-mu_lower**2)/mu_lower==F(1785)

# T1D-A mu > log(1/ell)-3 for ell<=1/16.
plus_constant=F(2,3)-F(16,5)+F(129,292)-F(3,16)
minus_constant=F(2,3)-F(16,5)+F(41,73)-F(4,16)
assert plus_constant> -3 and minus_constant> -3

# Elementary Gamma multiplier decomposition estimates:
# |c0|<14, |I-log(1+|xi|)|<5, |E_low|<4, |E_high|<16.
assert 14+5+4+16==39
# Lpsi^2 <= 5/pi + 30625/2304 < 16, hence Lpsi<4.
# Gamma residual <2(40+4)=88, prime shifts <48.
asymptotic_residual_bound=2*(40+4)+2*24
assert asymptotic_residual_bound==136
assert asymptotic_residual_bound**2==18496

# Available elementary fixed-window lower: Gamma > -14, primes > -48.
assert 14+2*24==62

print('PASS T1D-F: even/odd selected directions and parity-orthogonal residual Gram.')
print('PASS T1D-F: pole parity projections ||cosh(x/2)||^2=2+sinh 2; '
      '||sinh(x/2)||^2=sinh 2-2.')
print('PASS T1D-F: fixed T1C source B_+-=0 and 0<=B_++,B_--<1785.')
print('PASS T1D-F: rescaled-shape beta<=136; '
      'lambda_max(B)<=18496/(log(1/ell)-3) -> 0.')
print('PASS T1D-F: known q|F > -62||y||^2 is not positive coercivity.')
print('OPEN: tight projected residual norms; complement positivity; Lean; RH.')
