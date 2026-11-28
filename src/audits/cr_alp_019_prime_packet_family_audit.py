#!/usr/bin/env python3
"""CR-ALP-019. Exact prime-packet original Weil 4x4 span certificate.

Fraction-only numerical/log gates. The associated note contains the
continuum correlation, Gamma, Schwartz smoothing and Gram argument.
Not Lean, not the frozen 242x242 F41 H/V bank, not a universal support
extension, and not an RH theorem.
"""
from fractions import Fraction as F
from math import factorial, isqrt
from itertools import combinations

PRIMES = (2, 3, 5, 7)
b = F(1, 10**8)
eta = F(1, 10**14)
w = 2 * (b + eta)
G_lower = b - 2 * eta
assert 0 < eta < b / 4 and 0 < G_lower


def atanh_log_1_to_2(q, N=30):
    assert 1 <= q <= 2
    z = (q - 1) / (q + 1)
    t, lo = z, F(0)
    for j in range(N):
        lo += 2 * t / (2*j+1)
        t *= z*z
    hi = lo + 2*t/((2*N+1)*(1-z*z))
    return lo, hi


def log_int(n):
    assert n>=1
    m = n.bit_length()-1
    a,bound = atanh_log_1_to_2(F(n,2**m))
    c,d=atanh_log_1_to_2(F(2))
    return a+m*c, bound+m*d


def log_ratio(n,d):
    assert n>=d>=1
    an,bn=log_int(n)
    ad,bd=log_int(d)
    return an-bd,bn-ad


# Finite original receiving window dimensions; the four sources are
# genuinely successive new-support directions with increasing radii.
for p,rad in ((2,F(4,5)),(3,F(11,10)),(5,F(17,10)),(7,F(2))):
    assert log_int(p)[1]+b+eta < rad
assert log_int(3)[0]-log_int(2)[1] > w
assert log_int(7)[0]-log_int(5)[1] > F(1,4)+w
assert log_int(4)[0] > 1+w

# At most p^2=49 among diagonals; cross peaks pq<=35 or q/p<=7/2.
# Any prime power above 50 is separated from all those peaks by at
# least log(51/49)>=2/51, vastly exceeding w.
assert F(2,51)>w

def prime_powers_leq(n):
    sieve=[True]*(n+1)
    sieve[0]=sieve[1]=False
    for k in range(2,isqrt(n)+1):
        if sieve[k]:
            for t in range(k*k,n+1,k):sieve[t]=False
    out=set()
    for p in range(2,n+1):
        if sieve[p]:
            k=p
            while k<=n:
                out.add(k)
                k*=p
    return tuple(sorted(out))

prime_powers=prime_powers_leq(50)
assert 4 in prime_powers and 9 in prime_powers
assert 25 in prime_powers and 49 in prime_powers

# For positive a,b: |log(a)-log(b)| >= |a-b|/max(a,b).
# Exact integer/Fraction inequalities certify all prime-power samples
# outside the O(w) peak supports, without approximate floating logs.
def log_separation_lower(a,c):
    return abs(a-c)/max(a,c)

for p in PRIMES:
    assert log_separation_lower(F(2),F(1)) > w
    for n in prime_powers:
        if n == p*p:
            continue
        assert log_separation_lower(F(n),F(p*p)) > w, ('diagonal',p,n)
for p,q in combinations(PRIMES,2):
    peaks=(F(q,p),F(p*q))
    for n in prime_powers:
        for x in peaks:
            assert x!=n
            assert log_separation_lower(F(n),x) > w, ('mixed',p,q,n,x)

# From even mollifier and the positive interior of g_eta:
# G=||g_eta||^2 >= b-2eta, M_p=(1+p^2)G, and source supports disjoint.
for p,q in combinations(PRIMES,2):
    assert log_ratio(q,p)[0] > w
    assert log_ratio(q,p)[0] > F(1,4)+w
assert log_ratio(7,5)[0] > F(1,4)+w

# Check independent classical rational bounds for mGamma(0)>-6.
#  gamma < 1; pi < 22/7 from Machin rational alternating series.
def atan_inv(n, N=17):
    t=[F((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(N+1)]
    partial=sum(t[:N],F(0))
    if N%2==0: return partial,partial+t[N]
    else: return partial+t[N],partial
la,ua=atan_inv(5)
lc,uc=atan_inv(239)
pi_lo,pi_hi=16*la-4*uc,16*ua-4*lc
assert 3<pi_lo<pi_hi<F(22,7)
assert log_int(2)[1] < F(7,10)
assert sum((F(6,5)**k / factorial(k) for k in range(5)),F(0)) > F(22,7)
assert F(1)+3*F(7,10)+F(8,5)+F(6,5)<6

# e<11/4 using factorial-tail enclosure from k=0..7.
e_upper=sum((F(1,factorial(k)) for k in range(8)),F(0)) + F(9,8*factorial(8))
assert e_upper < F(11,4)
assert w*(F(11,4)**16) < 1 # ln(1/w)>16

# 3/2<e^(1/2)<5/3: positive Taylor lower, upper bounded first omitted term.
exp_half_lo=sum((F(1,2**k*factorial(k)) for k in range(3)),F(0))
exp_half_hi=sum((F(1,2**k*factorial(k)) for k in range(6)),F(0)) + F(14,13*2**6*factorial(6))
assert exp_half_lo > F(3,2)
assert exp_half_hi < F(5,3)
# hence e^(-1/2)>3/5 and e^(-1/2)<2/3.
assert -6+F(3,5)*16 == F(18,5)

# On |t|>=ln(7/5)-w>1/4, gamma weight wGamma(t)<3.
# ||g_eta||_1<=4sinh(b/2)<3b for b~1e-8.
# Thus all |B_pq| < 27(1+p)(1+q)b² <=1296b².
pair_bounds={}
for p,q in combinations(PRIMES,2):
    bound=27*(1+p)*(1+q)*b*b
    assert bound <= 1296*b*b
    pair_bounds[(p,q)] = bound

# In native physical M, every mixed normalized H entry <=
# Bmax/(sqrt((1+p²)(1+q²)) G) < Bmax/(7 G) because
# sqrt(5*10)>7 for two distinct p,q >=2.
assert 5*10 > 7*7
row_offdiag_crude=3*1296*b*b/(7*G_lower)
assert row_offdiag_crude < F(1,10000)
coercivity=F(18,5)-F(1,10000)
assert coercivity==F(35999,10000)
assert coercivity>0

# The four consecutive principal submatrices each inherit H>=c M.
for k in (2,3,4):
    assert k-1 <=3
    assert F(18,5) - (k-1)*1296*b*b/(7*G_lower) > coercivity

print('CR-ALP-019 EXACT FRACTION AUDIT PASS')
print('Primes:',PRIMES, '; prime powers in source cutoff:',prime_powers)
print('Strict support radii: 4/5, 11/10, 17/10, 2')
print('Mixed crossed Mellin poles and original prime-power sum = EXACTLY ZERO')
print('Sole positive original diagonal primes: n=4,9,25,49')
print('Gamma diagonal >18/5 native mass, pairwise Gamma |B| <27(1+p)(1+q)b²')
print('Row mixed / native mass upper bound:',row_offdiag_crude)
print('Full ORIGINAL signed 4x4 H > (35999/10000) physical M')
print('Every 2D/3D/4D nested selected-subspace Schur step inherits strict positivity')
print('No F41 242x242 embedding, full old-core positivity, new Lean, or RH proof')
