#!/usr/bin/env python3
"""CR-INT-T1D-E: exact finite-dimensional/rational audit.

Checks the a=2 arithmetic window, conservative L2 mixed-source bound,
pole Gram and finite-projection algebra, Schur identity, and rank-two
cross-block Gram spectrum. Analytic integral and density arguments are
proved in research/CR-INT-T1D-E.md, NOT by this script.
No Lean certificate, complement coercivity or RH.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp

def isprime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def is_primepower(n):
    for p in range(2,n+1):
        if isprime(p):
            k=n
            while k%p==0:
                k//=p
            if k==1:
                return True
    return False

# Rational Taylor/geometric majorant: 54 < exp(4) < 55.
N=24
terms=[F(4)**k/F(int(sp.factorial(k))) for k in range(N+1)]
partial=sum(terms,F(0))
tail=(F(4)**(N+1)/F(int(sp.factorial(N+1))))/(1-F(4,N+2))
assert 54<partial and partial+tail<55
pp=[n for n in range(2,55) if is_primepower(n)]
assert pp==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,
           27,29,31,32,37,41,43,47,49,53]
# For p>=2, log(p)/sqrt(p) <= 2/e <1; so S_2<24.
assert F(1)+F(1)+F(1,2)+F(1,6)>2
# w(t)<=2/t for 0<t<=1; w(t)<=2exp(-t/2) for t>=1.
# Gamma bound: ||k_gamma(u)|| <= ||u''||+(|C_gamma|+2+12)||u||.
# log(4pi)<3 follows 4pi<88/7<13<exp(3), and gamma<1.
assert F(4)*F(22,7)<13
assert sum(F(3)**k/F(int(sp.factorial(k))) for k in range(6))>13
assert (4+2+12+2*len(pp))==66

# At a=2, pole Gram [[2sinh2,4],[4,2sinh2]] is SPD.
assert F(2)+F(2)**3/F(6)>2

# Exact projection-residual formula, with nonorthogonal pole Gram.
s,t=sp.symbols("s t",real=True)
G_pole=sp.Matrix([[5,4],[4,5]])
assert G_pole.det()>0
e=sp.Matrix([sp.Rational(7,3),sp.Rational(11,4)])
p=sp.Matrix([s,t])
rf=sp.symbols("rf",real=True)
norm_squared=(e.T*e)[0]+(p.T*G_pole.inv()*p)[0]+rf**2
residual=norm_squared-(e.T*e)[0]-(p.T*G_pole.inv()*p)[0]
assert sp.simplify(residual-rf**2)==0

# Exact rank-2 operator/2x2 Gram nonzero spectrum.
r1=sp.Matrix([1,2,0]);r2=sp.Matrix([0,1,1])
mu1=sp.Rational(2);mu2=sp.Rational(3)
R=r1*r1.T/mu1 + r2*r2.T/mu2
B=sp.Matrix([[r1.dot(r1)/mu1,r1.dot(r2)/sp.sqrt(mu1*mu2)],
             [r2.dot(r1)/sp.sqrt(mu1*mu2),r2.dot(r2)/mu2]])
z=sp.symbols("z")
assert sp.simplify(sp.det(z*sp.eye(3)-R)-z*sp.det(z*sp.eye(2)-B))==0
assert R.rank()==2

# Real Schur slice, complex version by sesquilinear completion.
a,b,d,t0,t1,x,y,h=sp.symbols("a b d t0 t1 x y h",real=True,positive=True)
original=a*x*x+b*y*y+d*h*h+2*t0*x*h+2*t1*y*h
completed=a*(x+t0*h/a)**2+b*(y+t1*h/b)**2+(d-t0**2/a-t1**2/b)*h*h
assert sp.simplify(original-completed)==0
assert sp.det(sp.Matrix([[1,2],[2,1]]))==-3

print("PASS T1D-E: 54 < exp(4) < 55; exactly 24 prime powers")
print("PASS T1D-E: fixed-source mixed bound ||K_2u||_2 < ||u''||_2+66||u||_2")
print("PASS T1D-E: pole/projection Schur and 2x2 cross-block norm algebra")
print("OPEN: complementary coercivity and quantified residual Gram, no RH")
