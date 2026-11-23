#!/usr/bin/env python3
"""CR-ALP-014 exact source-scale and normalized-cycle gate regression.

Proof status: finite algebra checks only; no Lean compilation, no actual
D7 operators, no RH or original-Weil operator theorem.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as S

# Positive 3-word frame Gram can be realized by two disjoint 3D sectors,
# one ramified and one complementary, plus 3(d-2) orthogonal columns.
a,b,c,t=F(1,5),F(1,6),F(1,7),F(1,10)
C0=F(3)
assert a*b*c==F(1,210)

# These rational sigma values are illustrative *admissible error parameters*,
# not sigma=e^(-x^(eps/8)) arising at any named arithmetic input x.
for p,sigma in [(3,F(1)),(5,F(1)),(11,F(1,2)),
                (101,F(1,7)),(1009,F(1,22))]:
    d=(p-1)**2; r=(p-2)**2; q=2*p-3
    Gamma=F(p*p)
    Dstar=d*Gamma**2
    Dall=3*d
    de=F(q,d)
    assert d-r==q and 0<q<d and de<=sigma**2
    assert F(r,d)==F(p-2,p-1)**2
    assert d>=1/sigma**2
    assert F(d*Gamma**2,Dstar)==1  # ChProviso comparator in [1/C0,C0]
    assert Dall<=Dstar               # full target dimension LinkS condition
    beta=min(F(1),C0*Gamma**2/Dstar)
    assert beta<=C0**2*sigma**2
    # All three offdiagonal HS^2 are <=1/(20 Gamma^2).
    hs_max=(t*t + max(a*a,b*b,c*c))/Gamma**2
    assert hs_max==F(1,20)/Gamma**2
    assert hs_max<=F(d*d,Dstar)

    # Derived better C-normalized cycle-error upper bound, Eq (14.2).
    # E=6(1+sigma) sigma sqrt(beta+sigma sqrt(beta)).
    beta_sym=S.Rational(beta.numerator,beta.denominator)
    sig=S.Rational(sigma.numerator,sigma.denominator)
    E2=S.simplify(36*(1+sig)**2*sig**2*(beta_sym+sig*S.sqrt(beta_sym)))
    Enew2=S.simplify(144*(C0*C0+C0)*sig**4)
    assert S.simplify(Enew2-E2)>=0
    # The weaker order-one B error budget in CR-ALP-013 cannot be
    # shown to vanish from this expression in the actual gamma regime.
    B13_sq=S.simplify(Dstar*E2)
    assert S.simplify(B13_sq-36*Gamma**2*sig**2/C0)>=0
    # Even the improved B-scale upper budget grows with sample p.
    Bnew_sq=S.simplify(Dstar*Enew2)
    assert S.simplify(Bnew_sq-
                      144*(C0*C0+C0)/C0*Gamma**2*sig**2)>=0

    # Two positive Gram 3x3 sectors with offdiagonal scaled by 1/Gamma.
    g=S.Matrix([[1,S.Rational(a/Gamma),S.Rational(b/Gamma)],
                [S.Rational(a/Gamma),1,S.I*S.Rational(c/Gamma)],
                [S.Rational(b/Gamma),-S.I*S.Rational(c/Gamma),1]])
    f=S.Matrix([[1,S.Rational(t/Gamma),S.Rational(t/Gamma)],
                [S.Rational(t/Gamma),1,S.Rational(t/Gamma)],
                [S.Rational(t/Gamma),S.Rational(t/Gamma),1]])
    assert g[:2,:2].det()>0 and S.simplify(g.det())>0
    assert f[:2,:2].det()>0 and f.det()>0
    # For sign ±, full 3-cycle is (1/1000 ± i/210)/Gamma^3;
    # signed phase is independent of p, despite tiny C and B error.
    T_scale=F(Dstar,d**3)/Gamma**3
    re_part=T_scale*t**3
    im_part=T_scale*a*b*c
    assert im_part/re_part==F(100,21)
    assert im_part==F(1,210*d*d*Gamma)
    # B-scale absolute difference = sqrt(Dstar)*im_part
    assert isqrt(d)**2==d
    B_actual=Gamma*isqrt(d)*im_part
    assert B_actual==F(1,210*(p-1)**3)
    # Absolute C-error sits below the source-level coarse budget.
    assert S.Rational(im_part.numerator,im_part.denominator)**2<=E2
    print(f'PASS p={p} d={d} q={q} q/d={de} sigma={sigma} Gamma={Gamma}'
          f' Bactual={B_actual} Gamma*sigma={Gamma*sigma}')

print('CR-ALP-014 EXACT SCALE/PHASE REGRESSION PASS')
print('The full-source scaling proves Gamma*sigma -> infinity;')
print('the computed *upper budgets*, NOT the true errors, can diverge.')
print('No arithmetic D7 instance, Lean kernel theorem, Weil-form intertwiner or RH claim.')
