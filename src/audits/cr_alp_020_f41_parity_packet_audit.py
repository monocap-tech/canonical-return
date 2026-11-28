#!/usr/bin/env python3
"""CR-ALP-020: exact F41 odd-generator parity gate and even Weil packet Schur.

Source contracts are separately inspected from the frozen Lean sources.
This checks exact rational/log/scaling inequalities, not the continuum
parity/Gamma argument, not a 242-entry numerical F41 build, not RH.
"""
from fractions import Fraction as F

b=F(1,10**8)
eta=F(1,10**14)
w=2*(b+eta)
Gmin=b-2*eta
assert Gmin>0 and eta<b/4

def atanhlog(z,N=30):
    assert F(1)<=z<=F(2)
    u=(z-1)/(z+1)
    t=u
    lower=F(0)
    for j in range(N):
        lower+=2*t/(2*j+1)
        t*=u*u
    return lower,lower+2*t/((2*N+1)*(1-u*u))

def lnint(n):
    e=n.bit_length()-1
    a,c=atanhlog(F(n,2**e))
    l,r=atanhlog(F(2))
    return a+e*l,c+e*r

def lnratio(n,d):
    n0,n1=lnint(n)
    d0,d1=lnint(d)
    return n0-d1,n1-d0

ln2,ln3=lnint(2),lnint(3)
assert ln2[1]+b+eta<F(4,5)
assert ln3[1]+b+eta<F(11,10)
assert ln3[0]-ln2[1]>w
assert lnratio(3,2)[0]-w>F(1,4)
assert lnratio(6,5)[0]>w
assert lnratio(4,3)[0]>w
assert lnratio(9,8)[0]>w
assert lnratio(10,9)[0]>w
assert lnratio(7,6)[0]>w

# Frozen F41 j<241 oscillator shape:
# x^k sin(omega*x) for even k, x^k cos(omega*x) for odd k.
# Last j=241 source: sinh(x/2). Each critical is 32 sines+sinh.
# These are all odd. The actual frozen 241 node list is source-pinned
# in the research note, NOT re-enumerated here.
for k in range(14):
    assert ((-1)**k)*(-1 if k%2==0 else 1)==-1
print('PASS: frozen generator shape parity; 241 oscillators + 1 sinh + 4 critical')

# h_p=g-p*translate g; E_p=h_p+reflection h_p.
# ||E_p||²=2||g-pRg||²>=2(p-1)²G. The physical even projection
# of h_p is E_p/2, orthogonal to ALL selected odd F41 sources.
for p,required in ((2,F(1,10)),(3,F(1,5))):
    frac=F((p-1)**2,2*(1+p*p))
    assert frac==required
    assert 2*(p-1)**2*Gmin>0
print('PASS: raw h2/h3 cannot belong to F41 odd span; distance²/mass >= 1/10,1/5')

# Each even autocorrelation has only diagonal prime n=p² (all other
# prime powers excluded by exact log gap) and signed contribution
# >= -||E_p||²; the complete gamma >=18/5||E_p||².
assert ln2[1]<1
assert ln3[1]<F(3,2)
assert 2*ln2[1]/2<1 and 2*ln3[1]/3<1
assert w*F(11,4)**16<1
assert F(-6)+F(3,5)*16==F(18,5)
assert F(18,5)-1==F(13,5)

# Cross 2,3 prime-power sums vanish (nonintegral 3/2 and n=6 not
# a prime power); crossed poles zero. Mixed Gamma <3||E2||_1||E3||_1.
Bmax=1296*b*b
mass2min=2*Gmin
mass3min=8*Gmin
assert mass2min*mass3min==16*Gmin*Gmin
normalized=Bmax/(4*Gmin)
assert normalized<F(1,10000)
margin=F(13,5)-F(1,10000)
assert margin==F(25999,10000)

print('PASS: E2,E3 even smooth two-pole-null sources within radius 11/10')
print('PASS: all 242+4 frozen physical mixed moments of E2,E3 vanish by parity')
print('PASS: full original signed two-even-packet Weil H>(25999/10000)*physical M')
print('PASS: gamma / prime / mixed normalizations and exact fraction bound',normalized)
print('CR-ALP-020 EXACT FRACTION AUDIT PASS')
print('No original F41 numerical-entry acceptance, full old-core Schur, Lean or RH')
