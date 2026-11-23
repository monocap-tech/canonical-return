#!/usr/bin/env python3
"""CR-INT-T1D: exact-rational audit of profile-uniform smooth-core bounds.

Checks the algebra and all rational comparison constants in the companion
T1D note. The analytic inputs (original signed Weil formula, correlation
support, Gamma integral inequalities) are proved in the note; this script
is not a Lean proof, does not integrate arbitrary test functions, and does
not establish global Weil positivity or RH.

Requires Python 3 and SymPy. Run: python3 canonical_return_t1d_uniform_profile.py
"""
from fractions import Fraction as F
from math import factorial
import sympy as sp

S = sp.sqrt(3)
d = sp.symbols('d', positive=True, real=True)
k = sp.symbols('kappa29', positive=True, real=True)
R = sp.symbols('R29', real=True)
J = sp.symbols('J', positive=True, real=True)

# Crude rational enclosures needed for the uniform theorem.
assert F(5, 3)**2 < 3 < F(7, 4)**2
l2lower = 2*(F(1, 3)+F(1, 81)+F(1, 1215))
assert l2lower == F(842, 1215)
assert 5*l2lower > F(17, 5)  # log(32) > 17/5
# Euler's constant: gamma < H4 - log(4) - 1/10.
H4 = sum((F(1, i) for i in range(1, 5)), F(0))
gamma_upper = H4 - 2*l2lower - F(1, 10)
assert gamma_upper == F(2903, 4860) < F(3, 5)
# pi < 22/7 and e^(13/5) > the fifth degree Taylor polynomial > 88/7.
taylor_e = sum(((F(13, 5)**j)/factorial(j) for j in range(6)), F(0))
assert taylor_e > F(88, 7)
# Thus C_gamma=log(4pi)+gamma < 13/5+3/5=16/5.

# e^(1/16) <= 1/(1-1/16)=16/15; no undesired prime-power
# can enter the logarithmic collision windows about 3, 9, 27.
assert F(16, 15) < F(29, 27)
for n0, expected in [(3, [3]), (9, [9]), (27, [26, 27, 28])]:
    possible = [m for m in range(1, 40)
                if F(15*n0, 16) < m < F(16*n0, 15)]
    assert possible == expected, (n0, possible)
assert 26 == 2*13 and 28 == 4*7  # neither is a prime power
# First new prime-power threshold beyond the exact 3,9,27 offsets.
assert F(29, 27) < F(27, 25)  # log(29/27) < log(27/25)
assert F(29, 27) < F(9, 8)

# The exact T1B 4-bump to 2D Mellin-pole-null pencil.
X = sp.Matrix([[1,0], [-4/S,1], [1,-4/S], [0,1]])
L = sp.Matrix([[S**j for j in range(4)], [S**(-j) for j in range(4)]])
assert sp.simplify(L*X) == sp.zeros(2)
M = sp.simplify(X.T*X)
assert sp.simplify(M-sp.Matrix([[sp.Rational(22,3),-8*S/3],[-8*S/3,sp.Rational(22,3)]])) == sp.zeros(2)
assert M.det() == sp.Rational(292,9)

# T1B generalized coefficient laws, with a_j allowed to have any sign.
c1p=(-83+76*S)/146
c2p=(15+32*S)/73
c3p=-(33+12*S)/146
c1m=(83+76*S)/146
c2m=(15-32*S)/73
c3m=(33-12*S)/146
assert all(bool(x) for x in [0<c1p, c1p<sp.Rational(1,2),
                               0<c2p, c2p<1, -sp.Rational(1,2)<c3p, c3p<0,
                               0<c1m, c1m<sp.Rational(3,2),
                               -1<c2m, c2m<0, 0<c3m, c3m<sp.Rational(1,2)])

# For ell<=1/16: g0/J > log(2/ell)-C-ell,
# |Aj/J| < ell, log(32)>17/5, C<16/5, log(3)>1.
# Hence BOTH selected generalized Rayleigh values have a universal bound.
mu_plus_lower=F(1,5)+F(129,292)-F(3,16)
mu_minus_lower=F(1,5)+F(41,73)-F(4,16)
assert mu_plus_lower-F(9,20)==F(5,1168)>0
assert mu_minus_lower-F(1,2)==F(17,1460)>0

# The spectral splitting remains positive even if phi changes sign.
# Delta=(10sqrt3 log3 +83a1-64sqrt3 a2+33a3)/73.
# Using sqrt3>5/3, sqrt3<7/4, log3>1, |aj|<1/16:
split_lower = (F(50,3)-F(83+112+33,16))/73
assert split_lower == F(29,876)>0
# Smallest of the 1,i,-1 normalized phase gaps is
# ((11-4sqrt3)/22)*Delta > (2/11)*Delta.
phase_lower = F(2,11)*split_lower
assert phase_lower==F(29,4818)>0

# Newly activated prime 29 at log(29/27): only the 0-3 Toeplitz
# coefficient changes by t3 -> t3-kappa29*J*Rphi(log(29/27)).
# The diagonal of X.T*T*X remains unchanged, the mixed entry shifts.
extra = sp.zeros(4)
extra[0,3] = extra[3,0] = -k*J*R
H_extra = sp.simplify(X.T*extra*X)
assert sp.simplify(H_extra-sp.Matrix([[0,-k*J*R],[-k*J*R,0]]))==sp.zeros(2)
mp,mm=M[0,0]+M[0,1],M[0,0]-M[0,1]
shift_plus=sp.simplify((-k*R)/mp)
shift_minus=sp.simplify((k*R)/mm)
assert sp.simplify(shift_minus-shift_plus-sp.Rational(33,73)*k*R)==0

print('PASS — CR-INT-T1D exact-rational audit')
print('real C_c^infty bump, arbitrary sign, support width ell <= 1/16:')
print('  mu_plus >',mu_plus_lower,'> 9/20')
print('  mu_minus >',mu_minus_lower,'> 1/2')
print('  mu_minus - mu_plus >',split_lower)
print('  each adjacent 1,i,-1 phase gap >',phase_lower)
print('first extra prime-power interaction at ell = log(29/27)')
print('when phi is nonnegative and ell exceeds the threshold,')
print('  gap gain = (33/73)*(log(29)/sqrt(29))*R_phi(log(29/27)) >= 0')
print('No global spectral lower bound, no RH or Lean claim.')
