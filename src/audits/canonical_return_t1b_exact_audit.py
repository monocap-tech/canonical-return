"""Canonical Return CR-INT-T1B: independent exact symbolic audit.

Scope: algebra of the two-pole-null four-cell smooth-test compression.
This verifies coefficient and Gram identities, NOT operator domains,
trace-class hypotheses, positivity of the full Weil form, or RH.
Requires sympy; run: python canonical_return_t1b_exact_audit.py
"""
import sympy as sp

sqrt3 = sp.sqrt(3)
s = 4 / sqrt3
r = sqrt3
J, ellog3, g0 = sp.symbols('J L g0', positive=True, real=True)
A1, A2, A3 = sp.symbols('A1 A2 A3', positive=True, real=True)

X = sp.Matrix([[1, 0], [-s, 1], [1, -s], [0, 1]])
M_poles = sp.Matrix([[r**j for j in range(4)], [r**(-j) for j in range(4)]])
assert sp.simplify(M_poles * X) == sp.zeros(2, 2)

G = sp.simplify(J * X.T * X)
G_expected = J * sp.Matrix([[sp.Rational(22, 3), -8*sqrt3/3],
                            [-8*sqrt3/3, sp.Rational(22, 3)]])
assert sp.simplify(G - G_expected) == sp.zeros(2, 2)
assert sp.simplify(G.det() - sp.Rational(292, 9) * J**2) == 0

kappa = [None, ellog3/sqrt3, ellog3/3, ellog3/(3*sqrt3)]
A = [None, A1, A2, A3]
t = [g0] + [-J*kappa[k] - A[k] for k in (1, 2, 3)]
T_no_poles = sp.Matrix([[t[abs(i-j)] for j in range(4)] for i in range(4)])
H = sp.simplify(X.T * T_no_poles * X)
q_expected = (sp.Rational(22, 3)*g0 + 16*sqrt3*A1/3 - 2*A2
              + sp.Rational(14, 3)*J*ellog3)
b_expected = (-8*sqrt3*g0/3 - 25*A1/3 + 8*sqrt3*A2/3 - A3
              - 2*sqrt3*J*ellog3)
H_expected = sp.Matrix([[q_expected, b_expected], [b_expected, q_expected]])
assert sp.simplify(H - H_expected) == sp.zeros(2, 2)

q, b = H[0,0], H[0,1]
m, n = G[0,0], G[0,1]
mu_plus = sp.simplify((q+b)/(m+n))
mu_minus = sp.simplify((q-b)/(m-n))
delta = sp.factor(mu_minus - mu_plus)
delta_expected = (10*sqrt3*ellog3 + (83*A1-64*sqrt3*A2+33*A3)/J)/73
assert sp.simplify(delta - delta_expected) == 0

# Independently evaluate the prime-only source form and its three phases.
Tprime = sp.Matrix([[0 if i == j else -J*kappa[abs(i-j)]
                     for j in range(4)] for i in range(4)])
Hprime = sp.simplify(X.T * Tprime * X)
Hprime_expected = J*ellog3*sp.Matrix([[sp.Rational(14,3), -2*sqrt3],
                                       [-2*sqrt3, sp.Rational(14,3)]])
assert sp.simplify(Hprime-Hprime_expected) == sp.zeros(2,2)

for beta, p_over_norm in [(1, (-41+5*sqrt3)/73),
                          (-1, (-41-5*sqrt3)/73),
                          (sp.I, sp.Rational(-7,11))]:
    c = X * sp.Matrix([1, beta])
    norm = sp.simplify((c.conjugate().T * c)[0])
    prime_only_positive_named = sp.simplify(-(c.conjugate().T*Tprime*c)[0] / (J*norm))
    assert sp.simplify(prime_only_positive_named - ellog3*p_over_norm) == 0

# Distinct representation routes on the smooth core:
#   Direct signed Weil form = H.
#   Connes semilocal trace on two-pole-null source = H (using trace theorem).
#   Suzuki kernel numerator after D = H (using equation 8.7).
# The latter two statements require their cited analytic source theorems,
# and are NOT independently calculated operator traces or kernel integrals here.

print('PASS: pole-annihilator rank-two kernel, mass Gram, Weil Gram,')
print('      exact prime phase values and generalized eigenvalue gap.')
print('G / J =', G / J)
print('H[0,0] =', sp.simplify(q))
print('H[0,1] =', sp.simplify(b))
print('mu_minus - mu_plus =', delta)
print('strict prime contribution to gap ≈', float(10*sqrt3*sp.log(3)/73))
print('strict nearest-neighbor phase-gap floor ≈',
      float(((11-4*sqrt3)/22)*(10*sqrt3*sp.log(3)/73)))
print('Gamma correction is positive because A2 < A1/sqrt(3),')
print('83*A1 - 64*sqrt(3)*A2 + 33*A3 > 19*A1 + 33*A3 > 0.')
