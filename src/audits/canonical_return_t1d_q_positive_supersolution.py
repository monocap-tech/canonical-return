#!/usr/bin/env python3
"""CR-INT-T1D-Q — exact-rational positive supersolution for signed-Weil prime operator.

Audits original fixed-window 24 prime powers, ALL signed shift paths of
length <=3, exact interval endpoint order, and the pointwise inequality
A h < (203/25) h for the strictly positive bounded function
h = 1 + 20 A 1 + 8 A^2 1 (no continuum sampling).

Here A is the nonnegative-kernel SELF-ADJOINT auxiliary operator from
T1D-P, not a claim that the signed Weil prime quadratic form is positive.

Operator inference uses the standard weighted Schur test; prolate/Gamma
implications import prior source-level analytic lemmas (not Lean/CI).
Does NOT certify positivity of the omitted low-frequency Schur block,
Suzuki/Connes completed-domain equivalence, or RH.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import factorial, isqrt, exp

PP=[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]

def prime_base(n):
    for p in range(2,n+1):
        if n%p: continue
        if any(p%d==0 for d in range(2,isqrt(p)+1)): continue
        a=n
        while a%p==0:a//=p
        if a==1:return p
    return None
assert [n for n in range(2,55) if prime_base(n) is not None]==PP

def log_q(x,N=70):
    x=F(x); assert 1<=x<=2
    z=(x-1)/(x+1);t=z;a=F(0)
    for j in range(N):
        a+=2*t/(2*j+1);t*=z*z
    return a,a+2*t/((2*N+1)*(1-z*z))

def log_n(n):
    k=n.bit_length()-1
    a,b=log_q(F(n,2**k));c,d=log_q(2)
    return a+k*c,b+k*d

def sqrt_n(n,D):
    r=isqrt(n*D*D)
    return F(r,D),F(r+1,D)

D=10**12
Klo={};Kup={}
for n in PP:
    a,b=log_n(prime_base(n))
    s,t=sqrt_n(n,10**25)
    wlo,whi=a/t,b/s
    assert wlo>0
    Klo[n]=(wlo.numerator*D)//wlo.denominator
    Kup[n]=-((-whi.numerator*D)//whi.denominator)
    assert F(Klo[n],D)<=wlo<=whi<=F(Kup[n],D)

# Certified e^4 bracket, not rounded float comparisons.
S=sum((F(4)**k/factorial(k) for k in range(100)),F(0))
rem=(F(4)**100/factorial(100))/(1-F(4,101))
scale=10**24
Qlo=F((S.numerator*scale)//S.denominator,scale)
Qhi=F(-((-(S+rem).numerator*scale)//(S+rem).denominator),scale)
assert Qlo<S<S+rem<Qhi and 54<Qlo<Qhi<55

# W_k^- and W_k^+ step-function event tables, each scaled by D^k.
# For each positive path contribution, replace every positive coefficient
# kappa_n by its certified lower/upper integer grid representation.
events=[(defaultdict(int),defaultdict(int)) for _ in range(3)]
ACT=[(F(n),Klo[n],Kup[n]) for n in PP]+[(F(1,n),Klo[n],Kup[n]) for n in PP]
counts=[0,0,0]

def add(k,mn,mx,lo,hi):
    ratio=mx/mn
    if ratio>=Qhi:return False
    assert ratio<Qlo,('undecided path feasibility',k,ratio)
    l=(0,F(1)/mn);r=(1,F(1)/mx)
    elo,ehi=events[k-1]
    elo[l]+=lo;elo[r]-=lo
    ehi[l]+=hi;ehi[r]-=hi
    counts[k-1]+=1
    return True

for r1,l1,u1 in ACT:
    mn1=min(F(1),r1);mx1=max(F(1),r1)
    if not add(1,mn1,mx1,l1,u1):continue
    for r2,l2,u2 in ACT:
        rr2=r1*r2;mn2=min(mn1,rr2);mx2=max(mx1,rr2)
        if not add(2,mn2,mx2,l1*l2,u1*u2):continue
        for r3,l3,u3 in ACT:
            rr3=rr2*r3;mn3=min(mn2,rr3);mx3=max(mx2,rr3)
            add(3,mn3,mx3,l1*l2*l3,u1*u2*u3)
assert counts==[48,1290,28160],counts

allpoints=set()
for elo,ehi in events:
    allpoints.update(elo);allpoints.update(ehi)
assert len(allpoints)==6054,len(allpoints)

# Float is only a proposed permutation; the order is proved independently
# with certified disjoint rational bounds for EVERY adjacent pair.
def float_key(ev):
    kind,r=ev
    return float(r) if kind==0 else exp(4)*float(r)
points=sorted(allpoints,key=float_key)
def bracket(ev):
    kind,r=ev
    return (r,r) if kind==0 else (r*Qlo,r*Qhi)
for a,b in zip(points,points[1:]):
    a_lo,a_hi=bracket(a);b_lo,b_hi=bracket(b)
    assert a_hi<b_lo,('event order unresolved',a,b)
assert points[0]==(0,F(1)) and points[-1]==(1,F(1))

# h = 1+20 A1+8 A²1 >0 on I.  A h = A1+20 A²1+8 A³1.
# Bound positive numerator above and positive h below by separately
# summing upper/lower path-weight step functions; the inequality is
# checked exactly on every open event interval of the ORIGINAL source.
runlo=[0,0,0];runhi=[0,0,0]
lam=F(203,25)
min_slack=None; max_ratio_numeric=0;worst=None;min_h=None
for ev in points[:-1]:
    for j,(elo,ehi) in enumerate(events):
        runlo[j]+=elo[ev];runhi[j]+=ehi[ev]
        assert runlo[j]>=0 and runhi[j]>=runlo[j]
    l1,l2,_=runlo
    u1,u2,u3=runhi
    # Numerator common D³, h denominator D² (but scaled to D³).
    Ah_hi = u1*D**2+20*u2*D+8*u3
    h_lo = D**3+20*l1*D**2+8*l2*D
    slack=203*h_lo-25*Ah_hi
    assert slack>0,('supersolution failed',ev,F(Ah_hi,h_lo))
    if min_slack is None or slack<min_slack:
        min_slack=slack
    ratio=float(F(Ah_hi,h_lo))
    if ratio>max_ratio_numeric:
        max_ratio_numeric=ratio;worst=ev
    if min_h is None or h_lo<min_h:min_h=h_lo
# Tail events at the terminal boundary must cancel exactly.
for j,(elo,ehi) in enumerate(events):
    runlo[j]+=elo[points[-1]]
    runhi[j]+=ehi[points[-1]]
assert runlo==[0,0,0] and runhi==[0,0,0]
assert min_h>=D**3

# Weighted Schur test analytic inference: if a nonnegative symmetric
# kernel K satisfies Kh <=lambda h, then ||K||_2<=lambda. This also
# applies to countable sums of deterministic shifts. h bounded positive.
# The source Weil prime term is -A, hence >= -lambda * ||f||².
# Import (do not reprove) source Gamma/prolate bounds of T1D-M/O.
THETA=F(1,128)
B=2**21
q_margin=(1-THETA)*(21*F(693,1000)-F(49,10))-THETA*6-lam
schur=q_margin-F(155,182)
assert q_margin>F(1) and schur>F(1,2)
# B=2^20 does NOT close the same Schur sufficient estimate even
# at zero leakage using exact upper log2 bound and same lambda.
assert log_q(2)[1] < F(1733,2500)  # 0.6932
assert 20*F(1733,2500)-F(49,10)-lam < F(155,182)

# Imported prolate m_B = 4B/pi + O(log B) trace/defect estimates.
log2lo,log2hi=log_q(2)
assert log2lo>F(693,1000) and log2hi<F(347,500)
Du=F(2,9)*((21+2)*F(347,500)+F(3,2))
mlower=F(14*B,11)-Du/(1-THETA)
mupper=F(4*B,3)+Du/THETA
m_lo=mlower.numerator//mlower.denominator+1
m_hi=-(-mupper.numerator//mupper.denominator)-1
assert m_lo<m_hi and m_lo>0

print('PASS: exact 24 original prime powers, 48 signed shifts; 54<exp(4)<55')
print('PASS: exact path counts',counts,'and',len(points),'strict rational-verified endpoints')
print('PASS: weighted h=1+20 A1+8 A²1 strictly positive, exact event sweep Ah<(203/25)h')
print('PASS: ||A||_L2 <= 203/25 =',float(lam),'via analytic weighted Schur test')
print('Noncertifying diagnostic max upper/lower ratio:',format(max_ratio_numeric,'.12f'))
print('CERTIFIED conditioned source q(y,y) >',q_margin,'||y||² at B=2^21')
print('CONDITIONAL on T1D-I Schur ||R||<155/182: remainder >',schur,'||y||²')
print('Original prolate excluded rank:',m_lo,'<= m_B <=',m_hi)
print('METHOD LIMIT: B=2^20 not sufficient with fixed coarse Gamma and Schur constants.')
print('OPEN: positivity of excluded low-mode Schur block; all of F, domains, RH.')
