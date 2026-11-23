#!/usr/bin/env python3
"""CR-INT-T1D-J independently executable exact collision/rational audit.
Source analytic Gamma estimates are explained in research/CR-INT-T1D-J.md.
This checks a four-dimensional smooth-core compression, not full positivity.
"""
from fractions import Fraction as F
from math import isqrt, factorial
from pathlib import Path
import runpy
import sympy as S

# Independent T1C exact-fraction Gamma and selected-eigenvalue regression.
T=runpy.run_path(str(Path(__file__).with_name("canonical_return_t1c_rational_enclosures.py")))
assert T["g0_over_j"].lo>F(679,100)
assert T["mu_plus"].lo>7 and T["mu_minus"].lo>7
ell=F(1,10000)
r=S.sqrt(3); s=4/r
X=S.Matrix([[1,0],[-s,1],[1,-s],[0,1]])
L=S.Matrix([[r**j for j in range(4)],[r**(-j) for j in range(4)]])
assert S.simplify(L*X)==S.zeros(2,2)
cp=list(X*S.Matrix([1,1])); cm=list(X*S.Matrix([1,-1]))
w=[S.Integer(1),-s,S.Integer(1)]
Ny=2+s*s
Np=sum(x*x for x in cp); Nm=sum(x*x for x in cm)
assert S.simplify(Ny-S.Rational(22,3))==0
assert S.simplify(Np-(44-16*r)/3)==0
assert S.simplify(Nm-(44+16*r)/3)==0

# Center pair (twice coefficient of log(3), coefficient of log(2)).
E=[(2*j-3,0) for j in range(4)]
R=[(2*j-3,1) for j in range(3)]
Q=[(-a,-b) for a,b in R]
def primebase(n):
    for p in range(2,n+1):
        if n%p: continue
        if any(p%q==0 for q in range(2,isqrt(p)+1)): continue
        v=n
        while v%p==0: v//=p
        if v==1:return p
    return None
pp={n:primebase(n) for n in range(2,55) if primebase(n)}
assert len(pp)==24
maxexp=1/(1-ell)
def hits(A,B):
    out=[]
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            k=a[0]-b[0]; assert k%2==0
            z=F(3)**(k//2)*F(2)**(a[1]-b[1])
            m=max(z,1/z)
            for n in pp:
                ratio=max(m/n,F(n)/m)
                if ratio<maxexp:
                    assert ratio==1,(a,b,n,ratio)
                    out.append((i,j,n))
    return out
assert hits(E,R)==[(0,0,2),(1,1,2),(2,2,2)]
assert hits(E,Q)==[(1,2,2),(2,1,2),(3,0,2)]
assert hits(R,Q)==[(1,2,4),(2,1,4)]
assert sorted(hits(R,R))==[(0,1,3),(0,2,9),(1,0,3),(1,2,3),(2,0,9),(2,1,3)]
assert all(F(1)!=F(3)**((a[0]-b[0])//2)*F(2)**(a[1]-b[1])
           for a in E for b in R)
assert all(F(1)!=F(3)**((a[0]-b[0])//2)*F(2)**(a[1]-b[1])
           for a in R for b in Q)
l2,l3=S.symbols("l2 l3",positive=True)
kappa={n:(l2 if pp[n]==2 else l3 if pp[n]==3 else S.Symbol("log"+str(pp[n])))/S.sqrt(n)
       for n in pp}
modes=[(E,cp),(E,cm),(R,w),(Q,w)]
P=S.zeros(4)
for i,(A,c) in enumerate(modes):
    for j,(B,d) in enumerate(modes):
        P[i,j]=S.simplify(sum((-kappa[n]*c[a]*d[b] for a,b,n in hits(A,B)),S.Integer(0)))
Dp=2-2*s+s*s; Dm=2+2*s+s*s
assert S.simplify(P-P.T)==S.zeros(4)
assert S.simplify(P[0,1])==0
assert S.simplify(P[0,2]+l2*Dp/S.sqrt(2))==0
assert S.simplify(P[1,2]+l2*Dm/S.sqrt(2))==0
assert S.simplify(P[0,3]+l2*Dp/S.sqrt(2))==0
assert S.simplify(P[1,3]-l2*Dm/S.sqrt(2))==0
assert S.simplify(P[2,2]-S.Rational(22,3)*S.Rational(7,11)*l3)==0
assert S.simplify(P[3,3]-P[2,2])==0
assert S.simplify(P[2,3]-2*s*l2/2)==0

# Conservative exact rational consequences of full Gamma source bounds.
# w(t)<1 at t>=log 2; w(t)<2 at t>=log(7/5);
# w(t)<3 at t>=log(5/4). Support geometry is checked above.
assert 12*12<7*7*3                 # 2<s<7/3
assert S.Rational(7,4)**2>3 and S.Rational(5,3)**2<3
assert F(16,3)*F(22,3)>36
assert (F(44,3)+F(80,9))*F(22,3)>169
assert sum((F(7,10)**k/F(factorial(k)) for k in range(4)),F(0))>2
assert F(7,5)**2<2
assert F(49,24)**2*F(5,7)<4
assert F(50,9)**2/F(5)<9
assert maxexp<F(16,15) and maxexp<F(15,14)
gamma_diag=(F(4)*F(7,3)+2)*ell/7
gamma_RQ=3*ell*F(13,3)**2/7
gamma_plus=2*ell*F(14,3)*F(13,3)/6
gamma_minus=2*ell*F(26,3)*F(13,3)/13
assert gamma_diag<F(1,5000)
assert max(gamma_RQ,gamma_plus,gamma_minus)<F(1,1000)
diag_lo=F(679,100)+F(7,11)*F(109,100)-F(1,5000)
cross_hi=F(2,3)*F(7,20)+F(1,1000)
assert diag_lo>F(149,20) and cross_hi<F(1,4)
bplus_hi=F(1,2)*F(25,9)/6+gamma_plus
bminus_hi=F(1,2)*F(109,9)/13+gamma_minus
assert bplus_hi<F(1,4) and bminus_hi<F(1,2)
assert diag_lo-cross_hi>F(36,5)
assert 2*F(1,2)**2/7==F(1,14)
assert F(36,5)-F(1,14)>7
assert 7-F(3,4)>6
print("PASS CR-INT-T1D-J: prime-2/3/4/9 exact masks, physical metric,")
print("complement q|V>36/5; Schur(V)>7; full q|E+V>6.")
print("OPEN: infinite-dimensional complement, Suzuki/Connes transport, RH.")
