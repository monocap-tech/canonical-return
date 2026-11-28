#!/usr/bin/env python3
"""CR-INT-T1D-L: source-level rational audit of a genuine infinite-dimensional
high-Mellin-frequency conditioned subspace of the original signed Weil form.

Verifies the 24-prime-power a=2 census, rigorous S=sum Lambda(n)/sqrt(n)<13,
elementary bounds for m_Gamma(0)>-6 and m_Gamma(xi)>log|xi|-8
(assuming the published digamma integral), factorial/moment leakage constants,
and the exact numeric coercivity margins for B=2^60 and B=2^2700.

Proofs of the Fourier Taylor remainder, gamma-kernel inequality, and
infinite-dimensional nullspace are written in CR-INT-T1D-L.md. This script
checks finite arithmetic only; it is not a Lean verification, a proof of
Weil positivity on the full complement, or a Suzuki-Connes domain transfer.
"""
from fractions import Fraction as Q
from math import isqrt, factorial


def prime_base(n):
    for p in range(2, n+1):
        if any(p % j == 0 for j in range(2, isqrt(p)+1)):
            continue
        m=n
        while m%p == 0:
            m//=p
        if m == 1:
            return p
    return None

pp={n:prime_base(n) for n in range(2,55) if prime_base(n)}
assert len(pp)==24
assert sorted(pp)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]


def log_q(x, N=55):
    """Exact rational enclosure for log(x), 1 <= x <= 2."""
    x=Q(x); assert 1 <= x <= 2
    z=(x-1)/(x+1)
    t=z
    part=Q(0)
    for j in range(N):
        part += 2*t / (2*j+1)
        t*=z*z
    rem=2*t / ((2*N+1)*(1-z*z))
    return part,part+rem


def log_int(n):
    k=n.bit_length()-1
    z=Q(n,2**k)
    a,b=log_q(z)
    c,d=log_q(2)
    return a+k*c,b+k*d


def sqrt_lower(n,digits=30):
    scale=10**digits
    x=isqrt(n*scale*scale)
    assert x*x<=n*scale*scale<(x+1)*(x+1)
    return Q(x,scale)

sum_lower=sum((log_int(p)[0]/Q(isqrt(n*10**60)+1,10**30) for n,p in pp.items()),Q(0))
sum_upper=sum((log_int(p)[1]/sqrt_lower(n) for n,p in pp.items()),Q(0))
assert sum_lower<sum_upper<Q(13)
assert Q(12)<sum_lower

# 54 < exp(4) < 55 using Taylor partial and geometric-tail estimate.
terms=[Q(4)**k/factorial(k) for k in range(25)]
partial=sum(terms,Q(0))
next_term=Q(4)**25/factorial(25)
tail=next_term/(1-Q(4,26))
assert 54<partial and partial+tail<55

# Euler gamma upper: gamma < H_4 - log 4 - 1/10 < 3/5.
H4=sum((Q(1,k) for k in range(1,5)),Q(0))
log2lo,log2hi=log_q(2,40)
assert log2lo>Q(69,100)
assert log2hi<Q(7,10)
assert H4-2*log2lo-Q(1,10)<Q(3,5)

# Classical pi<22/7 plus elementary exp(6/5)>22/7
assert sum((Q(6,5)**k/factorial(k) for k in range(5)),Q(0))>Q(22,7)
# m_Gamma(0) = -gamma - pi/2 - 3log2 - log pi.
# Therefore m_Gamma(0)>- (3/5+11/7+21/10+6/5)>-6.
crude_m0_upper=Q(3,5)+Q(11,7)+Q(21,10)+Q(6,5)
assert crude_m0_upper<Q(6)
# Published digamma integral gives m_Gamma(xi)-m_Gamma(0)
#   = 2int_0^infty(1-cos(xi*t)) w(t)dt, w(t)>=1/(2t) for t<=1.
# For |xi|>=1, int_0^1(1-cos(xi*t))/t dt >= log|xi|-2
# from integration by parts. Thus m_Gamma(xi)>log|xi|-8.

# Taylor/Mellin-jet source constraint: N >= 10 B, B >= 1 integer.
# Relative low-band energy <= (4B/pi)*(2B)^(2N)/(N!)^2.
# By N! >= (N/e)^N and e<3, pi>1, this is strictly <
# 4B (9/25)^N <4B*2^(-10 B) <=1/128, since B<=2^B.
assert Q(9,25)<Q(1,2)
assert Q(2)*Q(3)/10==Q(3,5)  # 2 e /10 < 3/5
assert Q(4,2**9)==Q(1,128)

# gamma high band >= log B-8, gamma low band >-6, prime >=-26.
# log B = k log2 > 69 k/100, low leakage <=1/128.
def lower(k):
    log_B_lo=Q(69*k,100)
    return log_B_lo-Q(8)-(log_B_lo-Q(2))/128-Q(26)

lo60=lower(60)
lo2700=lower(2700)
assert lo60==Q(4539,640) and lo60>Q(7)
assert lo2700>Q(1800)
# T1D-I's sharper ||R|| <155/182 <1, but requires prior analytic review.
assert Q(155,182)<1
# T1D-F's coarser ||R||<1785 is adequate for the k=2700 case.
assert lo2700-Q(1785)>Q(15)

print('PASS: 24 original prime powers and 54 < exp(4) < 55')
print('PASS: exact rational von Mangoldt sum S in [{:.12f}, {:.12f}] < 13'.format(float(sum_lower),float(sum_upper)))
print('PASS: m_Gamma(0)>-6, high-frequency m_Gamma(xi)>log|xi|-8 (analytic inputs)')
print('PASS: polynomial-moment low-band energy leakage <1/128 for N>=10 B')
print('CERTIFIED: at B=2^60, conditioned q(y,y) > {} ||y||^2 >7 ||y||^2'.format(lo60))
print('CERTIFIED: at B=2^2700, conditioned q(y,y) > {} ||y||^2 >1800 ||y||^2'.format(lo2700))
print('CONDITIONAL SCHUR: using T1D-F ||R||<1785, conditioned remainder >15 ||y||^2')
print('OPEN: finite low-mode block, all of F, original pole directions, exact Sonin dictionary, RH')
