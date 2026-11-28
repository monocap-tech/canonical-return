#!/usr/bin/env python3
"""CR-ALP-015 exact rational-lattice phase-margin and certificate dominance.

This is a finite abstract source-schema model: it is NOT a D7 theta lift,
NOT new Lean kernel acceptance, NOT a signed Weil operator, NOT RH.
All assertions use exact integers/Fraction/SymPy, no numerical tolerances.
"""
from math import lcm
import sympy as S

I=S.I

def dot(u,v):
    return S.expand((u.conjugate().T*v)[0])

def gram(vs):
    return S.Matrix(3,3,lambda i,j:dot(vs[i],vs[j]))

def cycle(G):
    return S.expand(G[0,1]*G[1,2]*G[2,0])

v0=S.Matrix([1,0,0]); v1=S.Matrix([S.Rational(5,13),S.Rational(12,13),0])
assert dot(v0,v0)==1 and dot(v1,v1)==1

# A 3-word rational Q(i) lattice pairing family with a limiting phase,
# and all full diagonal defects exactly zero.
for n in (4,5,10,101,1000,10000):
    s=S.Rational(2*n,n*n+1)
    h=S.Rational(n*n-1,n*n+1)
    assert s*s+h*h==1
    grams=[]
    for sign in (-1,1):
        v2=S.Matrix([S.Rational(7,25),sign*I*S.Rational(24,25)*s,
                     S.Rational(24,25)*h])
        vs=[v0,v1,v2]
        assert all(dot(v,v)==1 for v in vs)
        G=gram(vs)
        assert G==G.conjugate().T
        for x in G:
            assert S.im(x).is_Rational and S.re(x).is_Rational
        assert G[:2,:2].det()>0
        expect_det=(S.Rational(12,13)*S.Rational(24,25)*h)**2
        assert S.simplify(G.det()-expect_det)==0
        assert G.det()>0
        assert G[0,1]==S.Rational(5,13)
        assert G[0,2]==S.Rational(7,25)
        assert G[1,2]==S.Rational(7,65)+sign*I*S.Rational(288,325)*s
        for i in range(3):
            for j in range(i+1,3):
                assert S.expand(G[i,j]*S.conjugate(G[i,j]))<S.Rational(1,3)
        cyc=cycle(G)
        assert cyc==S.Rational(49,4225)+sign*I*S.Rational(4032*n,21125*(n*n+1))
        # Exact algebraic bilinear pairing. With U=V=I and a_u=1,
        # Bu_u(phi,F)=F^T v_u phi=besS_u(phi,F)=emb_Q(i)(Bk_u(phi,F)).
        phi=1+2*I
        fieldvec=S.Matrix([2-I,-1+3*I,4])
        for v in vs:
            Bu=S.expand((fieldvec.T*v)[0]*phi)
            Bk=S.expand(sum(fieldvec[k]*v[k] for k in range(3))*phi)
            assert Bu==Bk and S.re(Bu).is_Rational and S.im(Bu).is_Rational
            assert dot(v,v)==1  # no source column (Bk) vanishes
        grams.append(G)
    assert grams[1] == grams[0].conjugate()
    lam=S.symbols('lam')
    assert S.expand(grams[0].charpoly(lam).as_expr()-grams[1].charpoly(lam).as_expr())==0
    # Generic InputS-like d_i=1,Dstar=3,Gamma^2=3,C0=1.
    G=grams[1]
    denoms=[]
    for x in G:
        denoms += [int(S.denom(S.re(x))),int(S.denom(S.im(x)))]
    N=lcm(*denoms)
    z=cycle(G)
    assert S.im(z)>0 and S.im(z)*N**3>=1
    assert (S.im(z)*N**3).is_Integer
    print(f'PASS Q(i) phase witness n={n}: Im cycle={S.im(z)}, '
          f'common entry denominator N={N}')

# Real positive, real negative, and zero cycles all remain possible
# under the same generic finite-frame/rationality/HS-bound schema.
v2positive=S.Matrix([S.Rational(7,25),0,S.Rational(24,25)])
v2negative=S.Matrix([S.Rational(7,25),-S.Rational(72,125),S.Rational(96,125)])
for v2,sgn in ((v2positive,1),(v2negative,-1)):
    G=gram([v0,v1,v2])
    assert all(G[i,i]==1 for i in range(3))
    assert G.det()>0 and G[:2,:2].det()>0
    assert all(G[i,j]**2<S.Rational(1,3) for i in range(3) for j in range(i+1,3))
    assert S.im(cycle(G))==0 and S.sign(cycle(G))==sgn
Gzero=gram([S.Matrix([1,0,0]),S.Matrix([0,1,0]),S.Matrix([0,0,1])])
assert Gzero==S.eye(3) and cycle(Gzero)==0

# The sigma-only COMPLETION ERROR UPPER-BOUND BUDGETS are too large
# to certify any phase when native Gamma grows. No lower bound on the
# unknown true error is being asserted.
for C0 in (1,3,5):
    for sn,sd in ((1,3),(1,10),(1,100)):
        sig=S.Rational(sn,sd)
        for Gamma in (100,1000,100000):
            d=sd**2//(sn**2)  # d >= sigma^-2, an integer
            Dstar=d*Gamma*Gamma
            beta=min(S.Integer(1),S.Rational(C0*Gamma**2,Dstar))
            assert S.Rational(d*Gamma**2,Dstar)==1
            assert beta<=S.Rational(C0*C0,d)
            E13sq=36*(1+sig)**2*sig**2*(beta+sig*S.sqrt(beta))
            Env2=(1+sig)**3/Dstar
            LB=36*sig**2*Gamma**2/(C0*(1+sig))
            assert S.simplify(E13sq/Env2-LB)>=0
            if LB>1:
                assert S.simplify(E13sq-Env2)>0
            E14sq=144*(C0*C0+C0)*sig**4
            LB14sq=(S.Rational(144*(C0*C0+C0),C0))*Gamma**2*sig**2/(1+sig)**3
            assert S.simplify(E14sq/Env2-LB14sq)>=0
            if LB14sq>1:
                assert S.simplify(E14sq-Env2)>0
print('PASS: sigma-only certificate budgets exceed maximum cycle envelope when Gamma*sigma >>1')
print('CR-ALP-015 EXACT AUDIT PASS')
print('Q(i) lattice examples are ABSTRACT, not instantiated arithmetic D7 period operators.')
print('A nonzero margin does not rescue the present sigma-only phase certificate.')
