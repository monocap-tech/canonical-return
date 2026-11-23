#!/usr/bin/env python3
"""CR-ALP-010: exact symbolic algebra for bilinear Bessel transport and a
restricted-pairing/full-period compression countermodel.

NOT Lean, NOT an actual D7 instance, NOT an original signed Weil operator,
and NOT an RH theorem. Requires SymPy.
"""
import sympy as S

I = S.I
adj = lambda M: M.conjugate().T

def cycle(maps):
    a, b, c = maps
    return S.simplify(S.trace(adj(a) * b * adj(b) * c * adj(c) * a))

def normalized_cycle(maps):
    mass = S.prod(S.trace(adj(a) * a) for a in maps)
    return S.simplify(cycle(maps) / mass)

# Finite models of the bilinear Bessel coefficients F^T T_u phi.
# These source/target spaces are complex Hilbert spaces with orthonormal bases.
T = [
    S.Matrix([1, 0, 0]),
    S.Matrix([1, 1, 0]),
    S.Matrix([1, I, 1]),
]
V = S.Matrix([[0, I, 0], [1, 0, 0], [0, 0, 1]])
assert S.simplify(adj(V) * V) == S.eye(3)
a = [2 + I, 1 - 2*I, -1 + I]       # wordwise nonzero scalars
u = [S.Integer(1), I, -I]           # wordwise source unitaries

# Pairing equality Bu_u(phi,F)=a_u besS_u(U_u phi,V F).
# Coefficient matrices obey A_u=a_u V^T T_u U_u; one COMMON target V.
A = [a[j] * V.T * T[j] * u[j] for j in range(3)]
assert cycle(T) == 1 + I
assert normalized_cycle(T) == (1 + I)/6
factor = S.prod(S.expand_complex(S.conjugate(z) * z) for z in a)
assert S.simplify(cycle(A) - factor * cycle(T)) == 0
assert S.simplify(normalized_cycle(A) - normalized_cycle(T)) == 0

# Restricted Bessel data P_u(t) E are identical for all t, where E:C -> C^2
# is exactly isometric. Full K_uv=P_u^* P_v can nonetheless reverse phase.
t = S.symbols('t', real=True)
c = [S.Matrix([1,0]), S.Matrix([1,1]), S.Matrix([1,I])]
E = S.Matrix([1, 0])
P = [v.row_join(t * S.conjugate(v)) for v in c]
assert all(S.simplify(P[j]*E-c[j]) == S.zeros(2, 1) for j in range(3))

expected_raw = (1 + t*t)**3 + I*(1-t*t)*(1+t*t)**2
expected_normal = S.Rational(1,4) + I*(1-t*t)/(4*(1+t*t))
assert S.simplify(S.factor(S.expand_complex(cycle(P))) - expected_raw) == 0
assert S.simplify(normalized_cycle(P) - expected_normal) == 0
assert S.simplify(cycle(P).subs(t, 0)) == 1 + I
assert S.simplify(cycle(P).subs(t, 2)) == 125 - 75*I
assert S.simplify(normalized_cycle(P).subs(t, 0)) == S.Rational(1,4) + I/4
assert S.simplify(normalized_cycle(P).subs(t, 2)) == S.Rational(1,4) - 3*I/20

print('CR-ALP-010 EXACT SYMBOLIC PASS')
print('BesselId pairing component preserves normalized nonreal 3-cell cyclic trace')
print('Restricted Bessel data do not determine full period-kernel phase')
print('Scope excludes BuRational, actual arithmetic D7 instances, Lean and RH')
