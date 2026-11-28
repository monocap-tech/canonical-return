#!/usr/bin/env python3
"""CR-ALP-018 exact rational old/new original-Weil 2x2 Schur gate.

Analytic dependencies (not computer-proven): complete original Weil Gamma
identity, Cauchy-Schwarz autocorrelation bound, L1 convolution contraction,
nonnegative mollifier, and elementary support geometry.

This auditor certifies all displayed numerical rational/log inequalities
and the precise 2x2 Schur consequence. NOT an F41 242x242 H/V proof,
not an entire old-core positivity theorem, not Lean or RH.
"""
from fractions import Fraction as Q


def log_atanh(x,N=32):
    assert Q(1) <= x <= Q(2)
    z=(x-1)/(x+1)
    t=z
    l=Q(0)
    for j in range(N):
        l+=2*t/(2*j+1)
        t*=z*z
    return l, l+2*t/((2*N+1)*(1-z*z))


def log_positive_ratio(n,d=1):
    assert n>=d>=1
    def Lint(v):
        m=v.bit_length()-1
        a,b=log_atanh(Q(v,2**m))
        c,e=log_atanh(Q(2))
        return a+m*c,b+m*e
    a,b=Lint(n)
    c,d0=Lint(d)
    return a-d0,b-c


b=Q(1,10**6)
eta=Q(1,10**12)
w=2*(b+eta)
assert eta<b/4
ln2=log_positive_ratio(2)
ln3=log_positive_ratio(3)
ln32=log_positive_ratio(3,2)
ln6=log_positive_ratio(6)
ln65=log_positive_ratio(6,5)
ln43=log_positive_ratio(4,3)
ln9=log_positive_ratio(9)
ln7=log_positive_ratio(7)

# Exact nested source support after centering:
# a_2=ln2+b+eta < 4/5, a_3=ln3+b+eta < 11/10.
assert ln2[1]+b+eta < Q(4,5)
assert ln3[1]+b+eta < Q(11,10)
assert Q(4,5)<Q(11,10)

# Distinct source packet supports disjoint (h2, h3), cross k(0)=0.
# Their difference-correlation can occur only near ln(3/2), ln6.
assert ln32[0]-w > Q(2,5)           # wGamma(t) <2 on cross support
assert ln43[0]>w                    # ln(3/2)+w < ln2
assert ln65[0]>w                    # ln6-w > ln5
assert ln7[0]>Q(19,10)              # actual IV011 a_old+a_new=1.9 cutoff
assert ln6[1]+w < ln7[0]            # also no n=7 etc
assert ln32[1]+w<ln2[0]

# Individual i=i supports have nonzero prime-power k(ln p²)=-p*||g||².
# At this narrow packet width no other prime powers lie in sample bands.
assert ln2[0] > w                   # no positive log n in central kernel
assert ln43[0]>w                    # ln4 has no n=3,5 collision
assert ln65[0]>w                    # ln6 is not a prime power
assert log_positive_ratio(9,8)[0]>w
assert log_positive_ratio(10,9)[0]>w

# Exact gamma lower bound of CR-ALP-017:
# mGamma(0)>-6, exp(-1/2)>3/5, ln(1/w)>11.
# w<1/100000 so ln(1/w)>ln(100000)>11 using e<11/4.
assert w<Q(1,100000)
assert (Q(11,4)**11)<Q(100000)
gamma_floor = -6 + Q(3,5)*11
assert gamma_floor==Q(3,5)

# g_eta remains equal to m_e*e^(u/2), m_e>=1 on [eta,b-eta].
# Hence ||g_eta||² >= b-2eta. For sources h_p=g-pT_{2ln p}g,
# distinct packets have disjoint support and ||h_p||²=(1+p²)||g_eta||².
qmin=b-2*eta
assert qmin>b/2
m2min=Q(5)*qmin
m3min=Q(10)*qmin
D2min=gamma_floor*m2min
D3min=gamma_floor*m3min
assert D2min>Q(29,10)*b
assert D3min>Q(59,10)*b

# ||g_eta||_1 <=||g_b||_1 =4sinh(b/2)<3b. Distinct h_p packet
# amplitudes 1,p imply ||h2||_1 <9b, ||h3||_1<12b.
# wGamma(t)=e^-t/2/(1-e^-2t)<2 for t>=2/5, since exp(4/5)>2.
assert 1+Q(4,5)+Q(4,5)**2/2 > 2
Bmax=2*(9*b)*(12*b)
assert Bmax==216*b*b
# Matrix determinant D2*D3-|B|²>0, with explicit conservative floors.
Schur_det_lower=D2min*D3min-Bmax*Bmax
assert Schur_det_lower>0

# Generalized metric M=diag(5||g_eta||²,10||g_eta||²) (exact
# orthogonality). sqrt(50)>7, so |B|/sqrt(m2m3) < Bmax/(7*qmin).
ratio_upper=Bmax/(7*qmin)
assert ratio_upper<Q(1,10000)
assert gamma_floor-ratio_upper > Q(5999,10000)

print('CR-ALP-018 EXACT FRACTION TWO-SOURCE SCHUR AUDIT PASS')
print('Old h_{2}: support inside radius 4/5; new h_{3}: inside 11/10')
print('Both original crossed Mellin pole moments zero')
print('Only original diagonal prime powers n=4,9 survive; mixed prime sum zero')
print('Gamma diagonal > 3/5 native L2 mass; mixed |B| <',Bmax,'[rational]')
print('2x2 original signed Hermitian H is PD; H >= 5999/10000*M')
print('No entire old collar, F41 H/V, full complement Schur, or RH proof')
