#!/usr/bin/env python3
"""CR-INT-T1D-M — exact-rational source-prolate audit.

Checks finite arithmetic / constant inequalities for a source-compatible
logarithmic time-band concentration operator. The trace, trace-defect,
Gamma/Ci and smooth source-splitting proofs are in the research note;
these functional-analytic claims are NOT machine-proved by this script.
No source-level full complement coercivity, Sonin equivalence, Lean or RH.
"""
from fractions import Fraction as F
from math import isqrt, factorial


def is_prime(n):
    return n >= 2 and all(n % j for j in range(2, isqrt(n)+1))


def prime_base(n):
    for p in range(2,n+1):
        if not is_prime(p):
            continue
        m=n
        while m % p == 0:
            m //= p
        if m == 1:
            return p
    return None

PP={n:prime_base(n) for n in range(2,55) if prime_base(n) is not None}
assert sorted(PP) == [2,3,4,5,7,8,9,11,13,16,17,19,23,25,
                     27,29,31,32,37,41,43,47,49,53]


def log_q(x, N=65):
    """Exact rational atanh interval of log x for 1<=x<=2."""
    x=F(x)
    assert 1<=x<=2
    z=(x-1)/(x+1)
    term=z
    total=F(0)
    for j in range(N):
        total+=2*term/(2*j+1)
        term*=z*z
    return total,total+2*term/((2*N+1)*(1-z*z))


def log_int(n):
    k=n.bit_length()-1
    xlo,xhi=log_q(F(n,2**k))
    tlo,thi=log_q(2)
    return xlo+k*tlo,xhi+k*thi


def sqrt_interval(n, digits=30):
    M=10**digits
    lo=isqrt(n*M*M)
    assert lo*lo<=n*M*M<(lo+1)**2
    return F(lo,M),F(lo+1,M)

slo=shi=F(0)
for n,p in PP.items():
    lglo,lghi=log_int(p)
    rtlo,rthi=sqrt_interval(n)
    slo+=lglo/rthi
    shi+=lghi/rtlo
assert slo<shi<F(61,5)     # S_2 <12.2, NOT just <13.
assert slo>F(12)

# With a=2, n contributes only for n<exp(4), and 54<exp(4)<55.
terms=[F(4)**j/factorial(j) for j in range(25)]
part=sum(terms,F(0))
remainder=(F(4)**25/factorial(25))/(1-F(4,26))
assert 54<part<part+remainder<55

log2lo,log2hi=log_q(2)
assert F(693,1000)<log2lo<log2hi<F(347,500)
# Classical pi<22/7 and exp(6/5)>22/7 give log(pi)<6/5.
assert sum((F(6,5)**j/factorial(j) for j in range(5)),F(0))>F(22,7)
# m_gamma(0)+gamma = -pi/2-3log2-log pi.
# Bound |Ci(X)| <=2/X by one integration by parts and tail integration.
# For X>=64, 2/X <=1/32. Hence m_gamma(X)>log X-49/10.
gamma_constant_bound=F(11,7)+3*F(347,500)+F(6,5)+F(1,32)
assert gamma_constant_bound<F(49,10)
# The global bound m_gamma(0)>-6 is imported from T1D-L.
assert F(3,5)+F(11,7)+3*F(7,10)+F(6,5)<6

# Time-band K_B on I=(-2,2): positive contraction, trace 4B/pi.
# Trace defect D=tr(K-K^2) = (2/pi^2)[int_0^4 sin^2(Bt)/t dt
#                                  +4 int_4^inf sin^2(Bt)/t^2 dt]
# <=(2/pi^2)[log(4B)+3/2], valid B>=1/4.
# If m eigenvalues >theta, m <=tr(K)+D/theta.
K=44
B=2**K
theta=F(1,128)
trace_upper=F(4*B,3)  # pi>3
trace_defect_upper=F(2,9)*(F((K+2)*347,500)+F(3,2))
dimension_upper=trace_upper+trace_defect_upper/theta
assert trace_defect_upper/theta<F(951)
assert dimension_upper<trace_upper+951

# For y perpendicular to the high-concentration eigenmodes,
# Fourier low-band energy <= theta||y||^2. Using
# m_gamma high > log(B)-49/10, low >-6, prime >=-2 S >-122/5:
coercivity= (1-theta)*(K*F(693,1000)-F(49,10))\
            -theta*6-F(122,5)
assert coercivity==F(15123,16000)
assert coercivity>F(9,10)

# Conditional import from T1D-I: ||R||<155/182 (separate analytic claim).
previous_cross_upper=F(155,182)
schur_margin=coercivity-previous_cross_upper
assert schur_margin==F(136193,1456000)>F(9,100)

# Independent coarse (T1D-F) alternative requires a larger B but
# does NOT rely on the more elaborate T1D-I finite-width enclosure.
K_COARSE=2640  # maybe sufficient? check below
coarse=(1-theta)*(K_COARSE*F(693,1000)-F(49,10))-theta*6-F(122,5)
if coarse<=1785:
    K_COARSE=2700
    coarse=(1-theta)*(K_COARSE*F(693,1000)-F(49,10))-theta*6-F(122,5)
assert coarse>1785

print('PASS T1D-M: 24 true prime powers, S in ({:.12f}, {:.12f}) < 61/5'.format(float(slo),float(shi)))
print('PASS T1D-M: rational Gamma constant proves m_Gamma(xi)>log|xi|-49/10 for |xi|>=64 (analytic Ci input)')
print('PASS T1D-M: sinc time-band trace, trace-defect, spectral-count bounds (analytic inputs)')
print('CUT: B=2^{}={}, theta=1/128'.format(K,B))
print('Finite excluded-prolate rank m < 4 B/3 + 951; exactly dimension bound is',str(dimension_upper))
print('CERTIFIED arithmetic bound: conditioned q(y,y) > {} ||y||^2'.format(coercivity))
print('CONDITIONAL on T1D-I ||R||<155/182: Schur remainder > {} ||y||^2'.format(schur_margin))
print('Coarse T1D-F import alternative: k={} gives {} >1785'.format(K_COARSE,coarse))
print('OPEN: actual finite smooth low-mode Schur block, complete F, Sonin/R1/R2/H2prime, RH')
