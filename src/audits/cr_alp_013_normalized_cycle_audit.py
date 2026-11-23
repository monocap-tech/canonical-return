#!/usr/bin/env python3
"""CR-ALP-013: exact finite arithmetic of the normalized phase-error gate.

Source contracts: frozen Lean LinkS_at, ChProviso_at, RamRank_at,
RamShareF. This is an exact *abstract* frame-model test; it neither
constructs a D7 theta lift nor compiles Lean nor asserts RH.

All assertions are exact fractions or symbolic SymPy identities.
"""
from fractions import Fraction as F
import sympy as S

I = S.I
a,b,c,t = F(1,5),F(1,6),F(1,7),F(1,10)
Ssum = a*a + b*b + c*c
assert Ssum == F(3889,44100)

# Distinct nonzero restricted and complementary cycles. The two
# complementary Gram matrices are positive and are conjugate.
for sign in (-1,1):
    g=S.Matrix([[1,S.Rational(a.numerator,a.denominator), S.Rational(b.numerator,b.denominator)],
                [S.Rational(a.numerator,a.denominator),1,sign*I*S.Rational(c.numerator,c.denominator)],
                [S.Rational(b.numerator,b.denominator),-sign*I*S.Rational(c.numerator,c.denominator),1]])
    f=S.Matrix([[1,S.Rational(1,10),S.Rational(1,10)],
                [S.Rational(1,10),1,S.Rational(1,10)],
                [S.Rational(1,10),S.Rational(1,10),1]])
    assert g==g.conjugate().T
    assert g[:2,:2].det()>0 and g.det()==S.Rational(40211,44100)>0
    assert f.det()>0 and f[:2,:2].det()>0
    assert S.simplify(g.charpoly().as_expr()-(S.Symbol('lambda')-1)*((S.Symbol('lambda')-1)**2-S.Rational(3889,44100))) == 0
    assert g[0,1]*g[1,2]*g[2,0] == sign*I/S.Integer(210)

# A 3-word positive frame into C^(3*d) can be built by a Cholesky/Gram
# factorization of a direct sum F ⊕ G_s ⊕ I_(3*(d-2)), with source
# C^d for each word. Each source has a fixed ramified C^r and complement
# C^q, with one active column in each sector; all other columns are
# wordwise mutually orthogonal. Cross Gram blocks have two active
# diagonal entries: t and G_s[i,j]. Diagonal K_ii=I_d.
# p is an illustrative odd prime; f=2 (consistent with rank-law parity).
for p in (3,5,11,101,1009):
    d=(p-1)**2
    r=(p-2)**2
    q=2*p-3
    assert d-r==q and r>0 and q>0
    de=F(q,d)
    assert F(r,d)==F(p-2,p-1)**2
    assert de<=F(2,p-1)  # Bernoulli source ratio bound
    Dstar=Dall=3*d
    gamma_squared=3
    C0=1
    assert F(d*gamma_squared,Dstar)==1  # ChProviso_at numeric premise
    assert Dall==Dstar
    # All offdiagonal HS^2 <= 1/20, whereas LinkS RHS=d^2/Dstar=d/3.
    assert t*t+max(a*a,b*b,c*c)==F(1,20)<=F(d*d,Dstar)
    e0=F(0)
    eta0=F(0)
    # Source q/d relative deficit is de; the compression error in
    # each off-diagonal cross block is |G_s[i,j]|^2, independent of sign.
    beta=min(F(1),F(C0*gamma_squared,Dstar))
    assert beta==F(1,d)
    for edge in (a,b,c):
        delta_normalized_sq=edge*edge/F(d*d)
        # Derived sharper bound when e0=0: s_ij^2 <= 2 de beta.
        assert delta_normalized_sq <= 2*de*beta
    # The restricted and full cycles have unchanging phases as p grows,
    # but the Dstar/(d^3) normalized cycles vanish in magnitude.
    factor=F(Dstar,d**3)
    restricted = factor*t**3
    imaginary_magnitude = factor*a*b*c
    assert imaginary_magnitude>0 and restricted>0
    assert imaginary_magnitude/restricted == F(1000,210)
    print(f'PASS rank-law p={p} d={d} r={r} q={q}, q/d={de}; '
          f'cyclic imaginary/real={imaginary_magnitude/restricted}')

# Source-level universal normalized estimate:
# sigma=exp[-x^(eps/8)] (0<sigma<=1);
# e0,de<=sigma^2, eta0<=sigma; beta=min(1,C0 Gamma^2/Dstar).
# s_ij^2 <= 2(de+sqrt(de*e0)) (beta+sqrt(e0*beta))
#           <=4 sigma^2 (beta+sigma*sqrt(beta)).
# For diagonal-comparable dimensions, beta bounds 1/d_i above.
for sn,sd in ((1,1),(1,2),(1,5),(1,10),(1,100)):
    sigma=S.Rational(sn,sd)
    for bn,bd in ((1,1),(1,4),(1,100),(1,10000)):
        beta=S.Rational(bn,bd)
        de=e0=sigma**2
        expr=2*(de+S.sqrt(de*e0))*(beta+S.sqrt(e0*beta))
        benchmark=4*sigma**2*(beta+sigma*S.sqrt(beta))
        assert S.simplify(expr-benchmark)==0
        assert S.simplify(benchmark-8*sigma**2)<=0
        # Error bound for C=Dstar*T/(d_i*d_j*d_k).
        budget=6*(1+sigma)*sigma*S.sqrt(beta+sigma*S.sqrt(beta))
        assert S.simplify(budget-12*S.sqrt(2)*sigma)<=0
print('CR-ALP-013 EXACT SYMBOLIC + RATIONAL AUDIT PASS')
print('Normalized cyclic values converge together; their phases need NOT agree.')
print('No evaluated arithmetic 3-word period operators, Lean compile, or RH claim.')
