#!/usr/bin/env python3
"""CR-ALP-012: exact rational source-rank and mixed-complement regression.

This verifies finite identities and samples. Not Lean, a D7 instantiation, RH,
or an arithmetic-to-original-Weil operator theorem. Requires SymPy.
"""
from fractions import Fraction as F
import sympy as S

# Frozen Lean RamRank_at: d_ram/d_full = (1-1/(T-1))**f.
# These are compatible exact integer samples, not actual theorem instances.
for T, f in [(5, 2), (101, 2), (11, 1)]:
    d = (T-1)**f
    ram = (T-2)**f
    q = d-ram
    ratio = F(q,d)
    assert F(ram,d) == F(T-2,T-1)**f
    assert 0 < q < d and 0 < ratio <= F(f,T-1)
    print('PASS rank sample:', 'T=',T,'f=',f,'full=',d,
          'ram=',ram,'complement=',q,'relative deficit=',ratio)

I=S.I
Fmat=S.Matrix([[1,S.Rational(1,10),S.Rational(1,10)],
               [S.Rational(1,10),1,S.Rational(1,10)],
               [S.Rational(1,10),S.Rational(1,10),1]])
a,b,c=S.Rational(1,5),S.Rational(1,6),S.Rational(1,7)
R=S.diag(1,0)
Q=S.eye(2)-R
norm2=lambda X: S.simplify(S.trace(X.conjugate().T*X))

for sign in (-1,1):
    G=S.Matrix([[1,a,b],[a,1,sign*I*c],[b,-sign*I*c,1]])
    K=lambda u,v:S.diag(Fmat[u,v],G[u,v])
    for u in range(3):
        assert K(u,u)==S.eye(2)
        assert norm2(K(u,u)-S.eye(2))==0
        assert norm2(Q*K(u,u)*Q)==1  # complements survive full diagonal isometry
    for u in range(3):
        for v in range(3):
            K0=R*K(u,v)*R
            omitted=K(u,v)-K0
            assert norm2(omitted)==norm2(K(u,v))-norm2(K0)
            if u!=v:
                # Here q_u=q_v=1, e_u=e_v=0, so (12.5) is delta^2<=2.
                assert norm2(omitted)<=2
                # Sample RamShareF de=1/2, d_u=d_v=2.
                assert norm2(omitted)<=F(1,2)*(2+2)
    fullcyc=S.simplify(S.trace(K(0,1)*K(1,2)*K(2,0)))
    restrictedcyc=S.simplify(S.trace(R*K(0,1)*R*K(1,2)*R*K(2,0)*R))
    assert fullcyc==S.Rational(1,1000)+sign*I/S.Integer(210)
    assert restrictedcyc==S.Rational(1,1000)
    print('PASS finite full near-isometry sign',sign,'cycle',fullcyc,
          'restricted',restrictedcyc)

print('CR-ALP-012 EXACT RANK + MIXED COMPLEMENT AUDIT PASS')
print('No new Lean, D7 instance, global signed Weil or RH certification')
