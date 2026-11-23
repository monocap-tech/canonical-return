#!/usr/bin/env python3
"""CR-ALP-003 finite table and smoothed-source numerical regression.

Validates selected identities, not the analytic mollifier convergence proof,
source theorem domains, Alpoge's estimate M, or Lean acceptance.
"""
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 45

def base_ramp(y, X):
    return max(mp.mpf('0'), X*mp.exp(-abs(y)/2)-mp.exp(abs(y)/2))

def weight_ramp(y, a, w):
    if a <= w:
        return mp.mpf('0')
    return w*base_ramp(y, a/w)

# For any finite cyclic group, each column belongs to |T2| translates.
for n in (2, 4, 8, 11):
    H = list(range(n))
    T2 = {z for z in H if z % 3 != 2}
    delta = Fraction(len(T2), n)
    for col in H:
        count = sum(1 for u in H if any((u+z)%n == col for z in T2))
        assert Fraction(count, n) == delta
    for qsign in (-1, 1):
        for dsign in (-1, 1):
            row = Fraction(1+qsign, 2)
            averaged = delta * row * (1-dsign)
            character_combo = delta*Fraction(1+qsign-dsign-qsign*dsign,2)
            assert averaged == character_combo

# Each weighted table cell is a single-source autocorrelation.
for aa, ww in [(mp.mpf(80),mp.mpf(4)), (mp.mpf(150),mp.mpf(15)),
               (mp.mpf(20),mp.mpf(25)), (mp.mpf(8),mp.mpf('0.25'))]:
    for n in (2, 3, 5, 9, 19, 37, 101):
        weighted = weight_ramp(mp.log(n), aa, ww)/mp.sqrt(n)
        expectation = max(mp.mpf(0), aa/mp.mpf(n)-ww)
        assert mp.almosteq(weighted, expectation,
                           rel_eps=mp.mpf('1e-38'),
                           abs_eps=mp.mpf('1e-38'))

# Even compact triangle mollifier used only for a numerical moment regression.
# The analytical proof in the note uses an arbitrary C-infinity mollifier.
X = mp.mpf(11)
L = mp.log(X)
for delta in [mp.mpf('0.18'), mp.mpf('0.06'), mp.mpf('0.02')]:
    rho = lambda v: (delta-abs(v))/(delta*delta) if abs(v)<=delta else mp.mpf('0')
    mom_plus=mp.quad(lambda v: rho(v)*mp.exp(v/2),[-delta,0,delta])
    mom_minus=mp.quad(lambda v: rho(v)*mp.exp(-v/2),[-delta,0,delta])
    assert mp.almosteq(mom_plus,mom_minus,rel_eps=mp.mpf('1e-37'))
    assert 1 < mom_plus < mp.exp(delta/2)
    pole=2*(X-1)*L*mom_plus*mom_minus
    assert pole > 2*(X-1)*L

print('CR-ALP-003 FINITE TABLE/RAMP/MOMENT REGRESSIONS PASS')
print('Scope excludes M, smoothing topology proof, Lean and operator-domain transfer')
