#!/usr/bin/env python3
"""CR-ALP-006: pure-exact-rational interval proof for two physical cells.

Only classical bounds 223/71 < pi < 22/7 and the elementary Euler
constant bounds H_N-log N-1/(2N)<gamma<H_N-log N-1/(2(N+1))
are external analytic premises. All arithmetic is Fraction; no floats
are used for any assertion. No Lean, theta-lift, RH, or operator claim.
"""
from fractions import Fraction as Q


def prime_power_base(n):
    for p in range(2, n + 1):
        if n % p:
            continue
        if any(p % d == 0 for d in range(2, p)):
            continue
        z = p
        while z < n:
            z *= p
        return p if z == n else None
    return None


def log_atanh_bounds(q, N=18):
    """Return exact enclosing rationals for log(q), 1<=q<=2."""
    assert Q(1) <= q <= Q(2)
    z = (q - 1) / (q + 1)
    total, term = Q(0), z
    for j in range(N):
        total += 2 * term / (2*j + 1)
        term *= z*z
    rem = 2 * term / ((2*N+1)*(1-z*z))
    return total, total+rem


def log_integer_bounds(n, N=18):
    assert n >= 1
    k = n.bit_length() - 1
    q = Q(n, 2**k)
    a,b = log_atanh_bounds(q,N)
    c,d = log_atanh_bounds(Q(2),N)
    return a+k*c, b+k*d


def log_ratio_bounds(a,b=1,N=18):
    la,ua = log_integer_bounds(a,N)
    lb,ub = log_integer_bounds(b,N)
    return la-ub,ua-lb


def expression_bounds(terms,N=18):
    low=high=Q(0)
    for c,a,b in terms:
        a0,a1=log_ratio_bounds(a,b,N)
        if c >= 0:
            low+=c*a0
            high+=c*a1
        else:
            low+=c*a1
            high+=c*a0
    return low,high


def euler_bounds(N=100,terms=18):
    H=sum((Q(1,j) for j in range(1,N+1)), Q(0))
    a,b=log_integer_bounds(N,terms)
    return H-b-Q(1,2*N), H-a-Q(1,2*(N+1))


# f(u)=exp(u/2) 1_[0,log9](u); g(u)=f(u-log9).
# Both masses are 8 and their supports meet only at a measure-zero point.
mass=Q(8)
assert mass > 0
assert Q(3,2)*Q(8,9)==Q(4,3)

# Prime term for the symmetrized cross autocorrelation:
# 3(n-1)/n for prime powers 2<=n<=9; (81-n)/(3n) for 9<n<=81.
prime_coeff_cross={}
for n in range(2,82):
    p=prime_power_base(n)
    if p is None:
        continue
    c=Q(3*(n-1),n) if n<=9 else Q(81-n,3*n)
    prime_coeff_cross[p]=prime_coeff_cross.get(p,Q(0))+c
assert prime_coeff_cross[3]==Q(16,3)
assert prime_coeff_cross[2]==Q(533,64)
assert len(prime_coeff_cross)==22

# Pole term (80/3)log9; full Gamma is
# -3log5 -(27/2)log(50/41)+(1/3)log9+(1/6)log(82/81).
# Combined pole/Gamma coefficient of log9 is exactly 27.
terms_cross=[(Q(27),9,1)]
terms_cross += [(-c,p,1) for p,c in prime_coeff_cross.items()]
terms_cross += [(Q(-3),5,1), (Q(-27,2),50,41), (Q(1,6),82,81)]
lo_cross,hi_cross=expression_bounds(terms_cross)
assert Q(-132426,1000000) < lo_cross < hi_cross < Q(-132425,1000000)

# Diagonal original Weil form: 16log9-prime-8EulerGamma
# -8log(pi)+10log5-16log4.
prime_coeff_diag={}
for n in range(2,10):
    p=prime_power_base(n)
    if p is None:
        continue
    prime_coeff_diag[p]=prime_coeff_diag.get(p,Q(0))+2*(Q(9,n)-1)
assert prime_coeff_diag=={2:Q(39,4),3:Q(4),5:Q(8,5),7:Q(4,7)}
terms_diag=[(Q(16),9,1)]+[(-c,p,1) for p,c in prime_coeff_diag.items()]
terms_diag += [(Q(10),5,1),(Q(-16),4,1)]
lo_diag,hi_diag=expression_bounds(terms_diag)
lo_gamma,hi_gamma=euler_bounds()
lo_diag-=8*hi_gamma
hi_diag-=8*lo_gamma

# Classical 223/71 < pi < 22/7.
lo_pi,_=log_ratio_bounds(223,71)
_,hi_pi=log_ratio_bounds(22,7)
lo_diag-=8*hi_pi
hi_diag-=8*lo_pi
assert Q(450,1000) < lo_diag < hi_diag < Q(456,1000)

# M=8I, H=[[d,c],[c,d]], c<0, d+c>0.
assert lo_diag + lo_cross > Q(317,1000)
assert lo_diag - hi_cross > Q(582,1000)

for label,low,high in [
    ('cross H12',lo_cross,hi_cross),
    ('diagonal H11',lo_diag,hi_diag),
    ('minimum ordinary eigenvalue',lo_diag+lo_cross,hi_diag+hi_cross),
    ('minimum generalized eigenvalue',(lo_diag+lo_cross)/mass,(hi_diag+hi_cross)/mass),
]:
    print(f'{label}: [{float(low):.12f}, {float(high):.12f}] (display only)')
print('CR-ALP-006 PURE RATIONAL INTERVAL PASS')
print('Conclusion: physically orthogonal nonsmooth sources, strictly negative')
print('mixed signed Weil pairing, positive 2D restriction.')
print('Smooth disjoint-support approximants exist analytically;')
print('no numerical two-word arithmetic theta-lift instance is claimed.')
