#!/usr/bin/env python3
"""CR-ALP-011 exact symbolic, finite frame/ramified complement regression.

Not a D7 arithmetic instance or Lean certificate; all assertions symbolic.
Requires SymPy. The analytic source-contract reasoning is in the paired note.
"""
import sympy as S

I = S.I
t = S.Rational(1, 10)
a, b, c = S.Rational(1, 5), S.Rational(1, 6), S.Rational(1, 7)

F = S.Matrix([[1, t, t], [t, 1, t], [t, t, 1]])
assert S.factor(F.det()) == (1-t)**2 * (1+2*t)


def gs(sign):
    return S.Matrix([[1, a, b], [a, 1, sign*I*c],
                     [b, -sign*I*c, 1]])


def k(sign, i, j):
    G = gs(sign)
    return S.diag(F[i,j], G[i,j])


Dstar = S.Integer(6)
Gam = S.sqrt(3)
d = S.Integer(2)
assert S.simplify(d*Gam**2/Dstar) == 1

for sign in [-1,1]:
    G = gs(sign)
    assert S.simplify(G-G.conjugate().T) == S.zeros(3)
    assert G[:2,:2].det() > 0
    assert S.simplify(G.det()) == S.Rational(40211,44100)
    assert G.det() > 0

    # Full 6D source Gram of the B_i:C² -> C^6, with each B_i
    # carrying one ramified and one complementary vector.
    full = S.zeros(6)
    for i in range(3):
        for j in range(3):
            full[2*i:2*i+2, 2*j:2*j+2] = k(sign,i,j)
    lam = S.symbols('lam')
    assert S.factor(full.charpoly(lam).as_expr()) == S.factor(
        F.charpoly(lam).as_expr()*G.charpoly(lam).as_expr())
    assert all(k(sign,i,i)==S.eye(2) for i in range(3))
    assert all(S.simplify(k(sign,i,j)-k(sign,j,i).conjugate().T)
               == S.zeros(2) for i in range(3) for j in range(3))

    # Literal LinkS_at-shaped numerical constraints:
    # d_i=2, Dall=Dstar=6, Gamma=sqrt3, C0=1, eta0=eta1=e0=0.
    for i in range(3):
        for j in range(3):
            if i!=j:
                hs=S.trace(k(sign,i,j).conjugate().T*k(sign,i,j))
                assert hs == t*t + S.Abs(G[i,j])**2
                assert hs <= d*d/Dstar

    # Canonical rank-one ramified isometric inclusion E(z)=(z,0).
    E=S.Matrix([1,0])
    assert E.conjugate().T*E == S.eye(1)
    R=E*E.conjugate().T
    assert R == S.diag(1,0)
    for i in range(3):
        assert R*k(sign,i,i)*R == R
        # Whole diagonal near-isometry is exact, but its complement is NOT small.
        assert S.trace((S.eye(2)-R)*k(sign,i,i)*(S.eye(2)-R)) == 1

    # The whole RESTRICTED Bessel Gram is unchanged with sign.
    assert all((E.conjugate().T*k(sign,i,j)*E)[0] == F[i,j]
               for i in range(3) for j in range(3))
    cyc = S.simplify(S.trace(k(sign,0,1)*k(sign,1,2)*k(sign,2,0)))
    assert cyc == t**3+sign*I*a*b*c
    assert S.simplify(cyc/8) == S.Rational(1,8000)+sign*I/S.Integer(1680)
    assert F[0,1]*F[1,2]*F[2,0] == S.Rational(1,1000)
    delta = k(sign,0,1)-R*k(sign,0,1)*R
    assert S.trace(delta.conjugate().T*delta) == a*a
    print('PASS sign=',sign,'full cyclic=',cyc,
          'normalized=',cyc/8,'restricted=1/1000')

# Both signs: pairwise full HS^2, complete frame spectrum, and ranks agree.
assert all(S.trace(k(1,i,j).conjugate().T*k(1,i,j))
           == S.trace(k(-1,i,j).conjugate().T*k(-1,i,j))
           for i in range(3) for j in range(3))
assert S.factor(gs(1).charpoly().as_expr()
                -gs(-1).charpoly().as_expr()) == 0
print('CR-ALP-011 EXACT SYMBOLIC PASS')
print('Zero FULL diagonal defect, fixed nontrivial restricted Bessel Gram,')
print('identical full spectral moments and HS^2, but opposite FULL cyclic orientation.')
print('No arithmetic D7, Lean, signed Weil operator or RH transfer claimed.')
