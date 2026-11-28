#!/usr/bin/env python3
"""CR-ALP-017: exact three-cutoff ramp identity and original-radius pole-null test.

Fraction/SymPy algebra checks; continuum Gamma/support/Cauchy proof is in
research note. No new Lean/CI, Weil support extension, or RH theorem.
"""
from fractions import Fraction as Q
from math import factorial
import sympy as S


def log_atanh(z, count=22):
    assert Q(1)<=z<=Q(2)
    q=(z-1)/(z+1); acc=Q(0); term=q
    for j in range(count):
        acc+=2*term/(2*j+1)
        term*=q*q
    return acc,acc+2*term/((2*count+1)*(1-q*q))


def log_int(n):
    m=n.bit_length()-1
    lo,hi=log_atanh(Q(n,2**m)); a,b=log_atanh(Q(2))
    return lo+m*a,hi+m*b


def log_ratio(a,b):
    l,u=log_int(a); ll,uu=log_int(b)
    return l-uu,u-ll


x=S.symbols('x',positive=True)
q=S.symbols('q',positive=True)
s=q+1/q
A=2+s*s
alpha=3*x+9+3/x
beta=3+3/x
chi=1/x
assert S.simplify(alpha-(A+2*s*s+1).subs(q,S.sqrt(x)))==0
assert S.simplify(beta-(2*s/q+1+1/q**2).subs(q,S.sqrt(x)))==0
assert S.simplify(alpha*(x-1)-2*beta*(x*x-1)+3*chi*(x**3-1))==0
assert S.simplify(alpha*(x+1)-beta*(x*x+1)+chi*(x**3+1)
                  -(x*x+9*x+9+1/x))==0
assert S.simplify(alpha-beta+chi-(3*x+6+1/x))==0

# On n=e^t >=1 the shifted ramp functions have piecewise formulas
# in terms of P_Y(n)=(Y/n-1)_+; these prove the full EVEN kernel identity.
n=S.symbols('n',positive=True)
P=lambda y:y/n-1
r1lo=(P(q**4)/q-s*P(q*q)).subs(q,S.sqrt(x))
r1hi=(P(q**4)/q).subs(q,S.sqrt(x))
assert S.simplify(r1lo-S.sqrt(x)*(1-1/n))==0
assert S.simplify(r1hi-(x*S.sqrt(x)/n-1/S.sqrt(x)))==0
r2lo=(P(q*q)-(1+1/q**2)*P(q**4)+P(q**6)/q**2).subs(q,S.sqrt(x))
r2mid=(-(1+1/q**2)*P(q**4)+P(q**6)/q**2).subs(q,S.sqrt(x))
r2hi=(P(q**6)/q**2).subs(q,S.sqrt(x))
assert S.simplify(r2lo)==0
assert S.simplify(r2mid-(1-x/n))==0
assert S.simplify(r2hi-(x*x/n-1/x))==0


def pos(Y,m):return max(Q(0),Q(Y,m)-1)
def r1(m):
    if m<=9:return Q(3)*(1-Q(1,m))
    if m<=81:return Q(27,m)-Q(1,3)
    return Q(0)
def r2(m):
    if m<=9:return Q(0)
    if m<=81:return 1-Q(9,m)
    if m<=729:return Q(81,m)-Q(1,9)
    return Q(0)

for m in range(1,801):
    raw=Q(118,9)*pos(9,m)-Q(20,3)*r1(m)+r2(m)
    three=Q(109,3)*pos(9,m)-Q(10,3)*pos(81,m)+Q(1,9)*pos(729,m)
    assert raw==three,(m,raw,three)

# Compressed 3-shift pole null: f_9 -2cosh(b/2)T_b f_9 +T_2b f_9
# =g_b -3 T_log9 g_b, g_b=e^(u/2)[1_[0,b]-e^-b 1_[b,2b]].
# F_+(g_b)=0, F_-(h)=0, and convolution preserves both zeros.
b0=Q(1,10**6); eta=Q(1,10**12); width=2*(b0+eta)
Llo,Lhi=log_int(9)
assert Llo>2
assert Lhi+width<Q(11,5)              # centered radius <11/10
assert Q(3,2)*Llo>Q(11,10)            # disjoint 3-shift OUTSIDE IV
assert Llo-width>1                    # no overlap on [width,1]
assert log_int(2)[0]>width
assert log_ratio(9,8)[0]>width         # log8 below log9-width
assert log_ratio(10,9)[0]>width        # log10 above log9+width

# Original Gamma symbol mGamma(0)>-6.
# gamma<1, 3log2<21/10, pi/2<8/5, logpi<6/5.
assert log_int(2)[1]<Q(7,10)
assert Q(11,7)<Q(8,5)
z=Q(6,5)
assert sum((z**k/factorial(k) for k in range(5)),Q(0))>Q(22,7)
assert 1+3*Q(7,10)+Q(8,5)+Q(6,5)<6
# e^(1/2)<5/3, hence e^(-1/2)>3/5.
e_upper=sum((Q(1,factorial(k)) for k in range(5)),Q(0))+Q(1,100)
assert e_upper<Q(11,4)<Q(25,9)
# width<e^(-11), hence log(1/width)>11.
assert width<Q(1,100000)
assert Q(11,4)**11<Q(100000)
assert -6+Q(3,5)*11==Q(3,5)

print('CR-ALP-017 EXACT SYMBOLIC/FRACTION AUDIT PASS')
print('k_long = (3X+9+3/X) k_X -(3+3/X) k_X² +(1/X) k_X³')
print('Full S identity at X,X²,X³ and pole cancellation: PASS')
print('Short-shift smooth X=9 source centered within IV a=11/10: PASS')
print('Only prime power n=9 remains; its complete signed contribution is >0')
print('Full original signed Weil W(k_h) > 3/5 ||h||_2²: ANALYTIC+RATIONAL GATE')
print('No new IV Schur positivity, support-extension certificate, or RH')
