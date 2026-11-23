#!/usr/bin/env python3
"""CR-ALP-001: exact algebraic checks, not a Lean or analytic-domain proof."""
import sympy as s

L = s.symbols('L', positive=True, real=True)
tau = s.symbols('tau', real=True)
y = s.symbols('y', real=True, nonnegative=True)
t = s.symbols('t', real=True)
X = s.exp(L)
h = X*s.exp(-y/2) - s.exp(y/2)

assert s.simplify(s.diff(h,y,2)-h/4)==0
assert s.simplify(h.subs(y,L))==0
assert s.simplify(s.diff(h,y).subs(y,L)+s.exp(L/2))==0

tri = s.exp(L/2)*(L-y) + s.integrate(
    (X*s.exp(-t/2)-s.exp(t/2))*(t-y), (t,y,L)
)/4
assert s.simplify(s.expand(tri-h))==0

F=2*s.integrate(h*s.cos(tau*y), (y,0,L))
w=(X+1-2*s.exp(L/2)*s.cos(tau*L))/(tau**2+s.Rational(1,4))
assert s.simplify(F-w)==0

for n in [1,2,3,9,27,31,100]:
    xx=s.Rational(30)
    test=(xx/s.sqrt(n)-s.sqrt(n)) if xx>n else s.Integer(0)
    rhs=s.sqrt(n)*((xx/s.Integer(n)-1) if xx>n else s.Integer(0))
    assert s.simplify(test-rhs)==0

ss=s.symbols('s', real=True, positive=True)
n=s.symbols('n',positive=True,real=True)
u=s.symbols('u',positive=True,real=True)
M=s.integrate((u/n-1)*u**(-ss-1), (u,n,s.oo),conds='none')
assert s.simplify(M-n**(-ss)/(ss*(ss-1)))==0

print('CR-ALP-001 SYMBOLIC PASS: Mellin ramp; prime normalization; tent identity; Fourier factorization')
print('Scope: exact formal algebra only; no Weil trace or RH theorem certified')
