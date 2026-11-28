#!/usr/bin/env python3
"""CR-INT-T1D-O: rational audit of the compressed prime-shift norm.

Mathematically proved in research/CR-INT-T1D-O.md:
1. Disintegration of L2((-2,2)) into shift residue fibers of P_I(T_t+T_-t)P_I.
2. Fiber path spectrum 2 cos(pi/(m+1)); Perron vector comparison.
3. Gamma multiplier and prolate trace/leakage estimates from T1D-M/N.
This program independently checks *finite* prime-power census, exact path
lengths, Perron rational certificates, outward log/sqrt enclosures, and the
consequent B=2^29 lower/Schur margins and rank bounds.

Not an integral/operator proof by code, prolate eigenspectrum, low-mode
Schur matrix, full Weil positivity, Lean certificate, or RH.
"""
from fractions import Fraction as F
from math import factorial, isqrt


def prime_base(n):
    if n<2:return None
    for p in range(2,n+1):
        if n%p: continue
        if any(p%q==0 for q in range(2,isqrt(p)+1)):continue
        m=n
        while m%p==0:m//=p
        if m==1:return p
    return None

PP={n:prime_base(n) for n in range(2,55) if prime_base(n) is not None}
assert sorted(PP)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,
                    31,32,37,41,43,47,49,53]

# Full original prime-power window: 54 < e^4 < 55 by a positive
# exponential series with rational tail, no decimal exponentials.
sum_exp=sum((F(4)**k/factorial(k) for k in range(25)),F(0))
tail=(F(4)**25/factorial(25))/(1-F(4,26))
assert F(54)<sum_exp<sum_exp+tail<F(55)

# Determine the maximum *essential* shift-chain length.
# Number m is the unique integer with n^(m-1)<e^4<n^m.
# The shift graph along fibers is a path with at most m vertices.
LENGTHS={}
for n in PP:
    m=2
    while n**m<54: m+=1
    # Our strict e^4 interval is (54,55); detect crossing it exactly.
    while n**m<=54: m+=1
    assert n**(m-1)<=54 and n**m>=55, (n,m)
    LENGTHS[n]=m
assert {n for n in PP if LENGTHS[n]==6}=={2}
assert {n for n in PP if LENGTHS[n]==4}=={3}
assert {n for n in PP if LENGTHS[n]==3}=={4,5,7}
assert {n for n in PP if LENGTHS[n]==2}==set(PP)-{2,3,4,5,7}

# A_m path-adjacency Perron vector comparison: if v>0 and A_m v <= c v,
# then max eig(A_m) <= c by the weighted infinity operator norm.
VECTORS={2:([1,1],F(1)),
         3:([70,99,70],F(99,70)),
         4:([8,13,13,8],F(13,8)),
         6:([20,36,45,45,36,20],F(65,36))}
for m,(v,c) in VECTORS.items():
    assert len(v)==m and all(x>0 for x in v)
    for j in range(m):
        av=(v[j-1] if j>0 else 0)+(v[j+1] if j+1<m else 0)
        assert F(av,v[j])<=c,(m,j,av,c)

# Certified exact-rational brackets for all von Mangoldt coefficients:
# kappa_(p^k) = log(p)/sqrt(p^k).
def log_q(x,N=65):
    x=F(x); assert 1<=x<=2
    z=(x-1)/(x+1); t=z; total=F(0)
    for j in range(N):
        total+=2*t/(2*j+1)
        t*=z*z
    return total,total+2*t/((2*N+1)*(1-z*z))

def log_n(n):
    k=n.bit_length()-1
    a,b=log_q(F(n,2**k)); c,d=log_q(2)
    return a+k*c,b+k*d

def sqrt_n(n,digits=28):
    u=10**digits; v=isqrt(n*u*u)
    assert v*v<=n*u*u<(v+1)*(v+1)
    return F(v,u),F(v+1,u)

sum_lo=sum_hi=F(0)
old_lo=old_hi=F(0)
for n,p in PP.items():
    a,b=log_n(p);c,d=sqrt_n(n)
    klo,khi=a/d,b/c
    upper=VECTORS[LENGTHS[n]][1]
    sum_lo+=upper*klo
    sum_hi+=upper*khi
    old_lo+=2*klo
    old_hi+=2*khi
assert sum_lo<sum_hi<F(55,4),(float(sum_lo),float(sum_hi))
assert old_lo>F(24)
assert old_lo-sum_hi>F(10)  # actual geometry saves >10 in norm bounds

# Form-lower bound uses the exact smooth pole-null source and imports
# two analytic Gamma inequalities from Suzuki and T1D-M:
# m_Gamma(ξ)>log|ξ|-49/10, |ξ|>=64, and m_Gamma(ξ)>-6 always.
# The compact-window prime operator has |q_prime(f)|<= C ||f||²,
# with independently checked C<55/4 (instead of 2*S<122/5).
K=29
B=2**K
THETA=F(1,128)
log2_lo,log2_hi=log_q(2)
assert F(693,1000)<log2_lo<log2_hi<F(347,500)
assert B>64
coercivity=(1-THETA)*(K*F(693,1000)-F(49,10))-THETA*6-F(55,4)
assert coercivity==F(164019,128000)>1
schur_margin=coercivity-F(155,182)
assert schur_margin==F(5005729,11648000)>F(2,5)

# Trace defect/rank imported analytic M/N inequalities, with rational checks.
# D_B <(2/pi²)(log(4B)+3/2), and 3<pi<22/7.
D_upper=F(2,9)*((K+2)*F(347,500)+F(3,2))
assert D_upper/F(127,128)<6
assert D_upper/THETA<655
rank_lower=F(14*B,11)-D_upper/(1-THETA)
rank_upper=F(4*B,3)+D_upper/THETA
min_m=rank_lower.numerator//rank_lower.denominator+1
max_m=-(-rank_upper.numerator//rank_upper.denominator)-1
assert min_m==683290247
assert max_m==715828537
assert min_m-4==683290243
assert F(22390054965613,max_m)>F(31000)

def outward(lo,hi,digits=12):
    unit=10**digits
    floor=lo.numerator*unit//lo.denominator
    ceil=-((-hi.numerator*unit)//hi.denominator)
    def dec(v): return str(v//unit)+'.'+str(v%unit).zfill(digits)
    return '['+dec(floor)+', '+dec(ceil)+']'

# The old cut B=2^44 is NOT silently retained, and the new rank remains huge.
print('PASS: 24 true prime powers, exact shift-chain path lengths (6,4,3,2)')
print('PASS: rational positive Perron vectors give ||S_log(n)|| <= 65/36,13/8,99/70,1')
print('Prime norm outward interval (12 d.p.):', outward(sum_lo,sum_hi), '< 55/4')
print('Old whole-line bound 2S outward interval:',outward(old_lo,old_hi), '; improvement >10')
print('B=2^29:',B,' conditioned Weil lower >',coercivity)
print('Conditional on T1D-I: Schur remainder >',schur_margin)
print('Prolate excluded rank bracket:',min_m,'<= m_B <=',max_m)
print('Additional constraints in physical complement >=',min_m-4)
print('PASS: >31000-fold reduction vs N old minimum excluded-rank count')
print('OPEN: omitted low-mode matrix, full complement, Suzuki-Connes domain transfer, RH')
