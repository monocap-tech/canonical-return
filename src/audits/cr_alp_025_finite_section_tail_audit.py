#!/usr/bin/env python3
"""CR-ALP-025: exact arithmetic gates for an explicit original-Weil tail.

Audits rational cutoffs, the quantitative sine-band leakage ceiling, its
positive two-pole-null tail margin, and exact finite Schur linear algebra.
It does NOT mechanize the continuum Fourier, digamma, or operator-domain
arguments. No Lean/kernel certification or fixed-collar positivity claim.
"""
from fractions import Fraction as F

A = F(11, 10)
R = 10**14
N = 250_000_000_000_000


def atanh_log_1_2(z, terms=40):
    assert F(1) <= z <= F(2)
    u = (z - 1) / (z + 1)
    term, lower = u, F(0)
    for j in range(terms):
        lower += 2 * term / (2*j + 1)
        term *= u*u
    upper = lower + 2*term / ((2*terms + 1)*(1-u*u))
    return lower, upper


def log_int(n):
    assert n >= 1
    m = n.bit_length() - 1
    lo, hi = atanh_log_1_2(F(n, 2**m))
    ln2lo, ln2hi = atanh_log_1_2(F(2))
    return lo + m*ln2lo, hi + m*ln2hi


# Previous stage's log symbol lower constant -8 and prime mass <8 give
# W(q) >= E_log(q) -16 ||q||² on the exact two-pole-null space.
assert 14 * log_int(10)[0] > 32
assert F(3) * (N+1) / (2*A) > 2*R  # q_{N+1}>2R, pi>3
# If q_n >=2R, |(exp(irx),sine_n)| <= 8/(3 sqrt(a) q_n).
# The analytic sum_n>N 1/n² <= 1/N then implies leakage <=
#    256*a*R / (9*pi³*N) < 256*a*R / (243*N).
leak_ceiling = F(256)*A*R/(243*N)
assert leak_ceiling == F(2816, 6075)
assert leak_ceiling < F(1, 2)
mu = 32*(1-leak_ceiling) - 16
assert mu == F(7088, 6075)
assert mu > 1
print('PASS: physical radius 11/10; cutoff R=10^14; sine rank N=250000000000000')
print('PASS: tail band leakage <', leak_ceiling, '< 1/2')
print('PASS: original signed two-pole-null W(q) >=', mu, '* ||q||² on Q_N')


# Exact finite-dimensional toy confirms Schur residual identity and shows
# why H>=0 alone is not a certificate. All values below are rationals.
def dot(v, w):
    return sum((x*y for x,y in zip(v,w)), F(0))


def matvec(mat, v):
    return [dot(row, v) for row in mat]


def prod(A0, B0):
    return [[sum((A0[i][k]*B0[k][j] for k in range(len(B0))), F(0))
             for j in range(len(B0[0]))] for i in range(len(A0))]


def transpose(A0):
    return [list(row) for row in zip(*A0)]


pole1 = [F(x) for x in (1, 1, 0, 0, 0)]
pole2 = [F(x) for x in (0, 1, 1, 0, 0)]
# Pole Gram [[2,1],[1,2]] has inverse [[2,-1],[-1,2]]/3.
pole_G_inv = [[F(2,3), F(-1,3)], [F(-1,3), F(2,3)]]

def pole_null_project(v):
    mom = [dot(pole1,v), dot(pole2,v)]
    coeff = matvec(pole_G_inv,mom)
    return [v[k]-coeff[0]*pole1[k]-coeff[1]*pole2[k]
            for k in range(len(v))]


e1=[F(x) for x in (1,0,0,0,0)]
e2=[F(x) for x in (0,0,0,1,0)]
e3=[F(x) for x in (0,0,0,0,1)]
vectors=[pole_null_project(e) for e in (e1,e2,e3)]
assert all(dot(pole, v)==0 for pole in (pole1,pole2) for v in vectors)
M=[[dot(v,w) for w in vectors] for v in vectors]
assert M == [[F(1,3),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(1)]]
# Check rank-two downdate M=I-F^T G^{-1}F.
Fmom=[[dot(pole,e) for e in (e1,e2,e3)] for pole in (pole1,pole2)]
downdate=prod(prod(transpose(Fmom),pole_G_inv),Fmom)
assert all(M[i][j]==(int(i==j)-downdate[i][j]) for i in range(3) for j in range(3))

# Independent abstract 2x2 block: p-block H, tail block D, cross b.
# Its exact Schur complement H-b²/D determines the form sign.
H=F(1,4);D=F(2);b=F(1)
assert D>0 and H>0 and H-b*b/D==F(-1,4)
q=-b/D
assert H+2*b*q+D*q*q==F(-1,4)  # positive 1x1 H is insufficient
H2=F(4); D2=F(3); b2=F(1)
assert H2-b2*b2/D2>0
# Our sufficient criterion conservatively uses only D>=mu:
assert H2-b2*b2/mu>0
print('PASS: rank-two physical-metric downdate and exact Schur residual identities')
print('PASS: positive finite matrix need not imply full form positivity')
print('CR-ALP-025 EXACT RATIONAL AUDIT PASS; continuum argument remains analytic')
