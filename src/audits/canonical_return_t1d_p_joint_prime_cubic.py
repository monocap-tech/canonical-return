#!/usr/bin/env python3
"""CR-INT-T1D-P joint original prime operator: rational path-sweep audit.
Enumerates ALL 48^3 signed 3-step shifts and checks physical-window
admissibility after every step, with exact sorted endpoint comparisons.
A floating-key sort is verified against exact rational exp(4) brackets.
Analytic Schur operator and imported Gamma/prolate facts are not Lean proofs.
"""
from fractions import Fraction as F
from collections import defaultdict
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

assert [n for n in range(2,55) if prime_base(n)]==PP

def log_rat(x,N=70):
    x=F(x); assert 1<=x<=2
    z=(x-1)/(x+1); t=z; a=F(0)
    for j in range(N):
        a+=2*t/(2*j+1)
        t*=z*z
    return a,a+2*t/((2*N+1)*(1-z*z))

def log_n(n):
    k=n.bit_length()-1
    a,b=log_rat(F(n,2**k));c,d=log_rat(2)
    return a+k*c,b+k*d

D=10**12
Kup={}
for n in PP:
    lo,hi=log_n(prime_base(n))
    scale=10**25
    root=isqrt(n*scale*scale)
    sqrtlo=F(root,scale)
    sqrthi=F(root+1,scale)
    klo,khi=lo/sqrthi,hi/sqrtlo
    low=(klo.numerator*D)//klo.denominator
    up=-((-khi.numerator*D)//khi.denominator)
    assert F(low,D)<=klo<=khi<=F(up,D)
    Kup[n]=up

# Analytically justified Taylor bounds e^4 in a tiny rational interval.
S=sum((F(4)**j/factorial(j) for j in range(100)),F(0))
tail=(F(4)**100/factorial(100))/(1-F(4,101))
M=10**24
Qlo=F((S.numerator*M)//S.denominator,M)
Qhi=F(-((- (S+tail).numerator*M)//(S+tail).denominator),M)
assert Qlo<S<S+tail<Qhi and 54<Qlo<Qhi<55

ACT=[(F(n),Kup[n]) for n in PP]+[(F(1,n),Kup[n]) for n in PP]
events=defaultdict(int)
count=0
# For a 3-step path with cumulative products r_j, all intermediate
# windows 1 < z*r_j < e^4 give z in (1/min(r_j),e^4/max(r_j)).
for r1,w1 in ACT:
    mn1=min(F(1),r1);mx1=max(F(1),r1)
    for t2,w2 in ACT:
        r2=r1*t2
        mn2=min(mn1,r2);mx2=max(mx1,r2)
        if mx2/mn2>=Qhi:continue
        for t3,w3 in ACT:
            r3=r2*t3
            mn=min(mn2,r3);mx=max(mx2,r3)
            ratio=mx/mn
            if ratio>=Qhi:continue
            assert ratio<Qlo,("ambiguous path vs exp(4)",ratio)
            wt=w1*w2*w3
            events[(0,F(1)/mn)]+=wt
            events[(1,F(1)/mx)]-=wt
            count+=1

def approx(ev):
    k,r=ev
    return float(r) if k==0 else exp(4)*float(r)

points=sorted(events,key=approx)

def exact_bracket(ev):
    k,r=ev
    return (r,r) if k==0 else (Qlo*r,Qhi*r)

# Verify every proposed adjacent ordering with exact rational arithmetic.
for a,b in zip(points,points[1:]):
    assert exact_bracket(a)[1]<exact_bracket(b)[0],("misordered",a,b)

running=top=0
for ev in points:
    running+=events[ev]
    assert running>=0
    top=max(top,running)
assert running==0
assert count==28160 and len(points)==6054
assert top<800*D**3
assert 800*10**3<93**3  # ||A||^3 <= ||A^3 1||_inf <800

# Consequences depend on earlier analytic Gamma and prolate inequalities.
theta=F(1,128)
q=(1-theta)*(22*F(693,1000)-F(49,10))-6*theta-F(93,10)
assert q==F(58771,64000)>F(9,10)
schur=q-F(155,182)
assert schur==F(388161,5824000)>F(1,16)

# Imported prolate trace/defect bounds, only arithmetic re-evaluated.
B=2**22
l2,h2=log_rat(2)
assert l2>F(693,1000) and h2<F(347,500)
defect=F(2,9)*(24*F(347,500)+F(3,2))
lower=F(14*B,11)-defect/(1-theta)
upper=F(4*B,3)+defect/theta
ranklo=lower.numerator//lower.denominator+1
rankhi=-(-upper.numerator//upper.denominator)-1
assert (ranklo,rankhi)==(5338202,5592921)

print("PASS: all 24 prime powers, 48 signed shifts")
print("PASS: exact windowed triple-path sweep",count,"paths,",len(points),"breakpoints")
print("PASS: joint cubic row-sum <",800,"computed upper",float(F(top,D**3)))
print("CERTIFIED operator norm ||A||_2 < 93/10")
print("CONDITIONAL Gamma/prolate B=2^22 Weil lower",q)
print("CONDITIONAL T1D-I Schur remainder",schur)
print("Prolate excluded-rank bracket",ranklo,rankhi)
print("OPEN: low-mode Schur, full complement, RH")
