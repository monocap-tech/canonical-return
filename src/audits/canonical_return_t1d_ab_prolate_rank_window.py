#!/usr/bin/env python3
"""CR-INT-T1D-AB: exact-rational rank bounds for the tiny-leakage prolate gate.

Analytic proof (in companion note): trace/defect lower count, convolution-smoothed
Fourier window, finite trigonometric approximation to band restriction, and
approximation-number/min-max conversion to a rank upper bound. This program
checks ALL rational numeric inequalities, not those functional-analytic steps.
No prolate eigenvalue enumeration, low-mode Weil positivity, Lean or RH.
"""
from fractions import Fraction as F


def atanh_log(q, terms=36):
    q = F(q)
    assert q >= 1
    z = (q - 1)/(q + 1)
    z2 = z*z
    power, partial = z, F(0)
    for j in range(terms):
        partial += 2*power/F(2*j+1)
        power *= z2
    return partial, partial + 2*power/(F(2*terms+1)*(1-z2))


def atan(q, terms):
    q = F(q)
    assert 0 < q < 1
    partial = sum(((-1)**j * q**(2*j+1) / F(2*j+1)
                   for j in range(terms)), F(0))
    e = q**(2*terms+1)/F(2*terms+1)
    return (partial, partial+e) if terms % 2 == 0 else (partial-e, partial)

# The original signed Weil source does NOT enter this independent time-band
# rank computation. Physical I=(-2,2); unitary Fourier transform;
# T_B=P_I F^-1 1_[-B,B] F P_I, theta=2^-38.
B = 2**21
THETA = F(1,2**38)
assert B == 2097152
p5 = atan(F(1,5),28)
p239 = atan(F(1,239),18)
pi_lo, pi_hi = 16*p5[0]-4*p239[1], 16*p5[1]-4*p239[0]
assert F(314159,100000) < pi_lo < pi_hi < F(22,7)
ln2_lo, ln2_hi = atanh_log(2,36)
assert F(693,1000)<ln2_lo<ln2_hi<F(347,500)

# T1D-N's exact trace/defect estimate, with new sharp rational substitution.
# Trace = 4B/pi; D_B<=2/pi^2 (log(4B)+3/2).
trace_lo = F(4*B)/pi_hi
trace_hi = F(4*B)/pi_lo
D_hi = F(2)/pi_lo**2 * (F(23)*ln2_hi + F(3,2))
count_lo = trace_lo - D_hi/(1-THETA)
assert F(2670173)<count_lo<F(2670174)
rank_lower = 2670174  # m integer, > 2670173
print('PASS T1D-AB: rigorous trace/defect count m >=',rank_lower)
print('   rational lower diagnostic:',float(count_lo),'(rounded for display only)')

# Convolution-window Fourier rank upper certificate.
# Let delta=1/120, a=2+delta. χ=1_[-2-delta/2,2+delta/2]
# convolved with k=14 uniform densities each supported in
# [-delta/(2k),delta/(2k)]. It satisfies 0<=χ<=1, χ=1 on I,
# support [-a,a] and its Fourier transform obeys
# |χhat(t)| <= (2/|t|) min(1,C/|t|)^k, C=2k/delta.
# Band-limited synthesis S: L2[-B,B]->L2(I) approximated by the
# degree-|n|<=N Fourier series of χ S f on [-a,a]. Approx rank<=2N+1.
# If Ω=(π N/a)-B > 0, Parseval/Cauchy and integral-test imply
# ||S-S_N||^2 <= 4B/(π² (2k+1) Ω) * (C/Ω)^(2k).
# We prove Ω>9300 with rational π and use π>3 in the error bound.
delta = F(1,120)
a = F(2)+delta
k = 14
N = 1346600
C = F(2*k)/delta
assert a==F(241,120) and C==3360
Omega_lo = pi_lo*N/a - B
assert Omega_lo>9300
O = F(9300)
error2 = F(4*B,9*(2*k+1))*F(1,O)*(C/O)**(2*k)
assert error2 < THETA
rank_upper = 2*N+1
assert rank_upper == 2693201
print('PASS: χ convolution window exactly 1 on I; Fourier degree N=',N)
print('PASS: πN/a-B>9300; certified squared S-approximation norm <theta')
print('   bound/theta =',float(error2/THETA),'(rounded display only)')
print('CERTIFIED excluded prolate rank interval:',rank_lower,'<=m<=',rank_upper)

# Fourier truncation with even χ preserves reflection and decomposes into
# cos(0..N) and sin(1..N) on physical I. Thus rank upper separately N+1/N.
# T1D-N's reflected-trace estimate |Tr(R T_B)|<4/3, hence each parity
# trace > (4B/pi)/2 - 2/3. Defect of each parity <= full D_B.
parity_count_lo = trace_lo/2 - F(2,3) - D_hi/(1-THETA)
assert parity_count_lo>1335084
assert parity_count_lo<1335085
parity_rank_lower=1335085
assert N+1 == 1346601
print('CERTIFIED parity even:',parity_rank_lower,'<=m_even<=',N+1)
print('CERTIFIED parity odd: ',parity_rank_lower,'<=m_odd <=',N)

# Physical Z=(P+E+L)^perp has codim 8, split evenly (4 and 4).
# Rank of prolate moment map restricted to smooth Z is in [m-8,m].
# Thus the minimum number of additional smooth moment constraints is
# >=m-8; maximum m. Parity lower m_sigma-4.
assert rank_lower-8==2670166
assert parity_rank_lower-4==1335081
print('CERTIFIED at least',rank_lower-8,'independent new moment constraints in Z')
print('CERTIFIED at least',parity_rank_lower-4,'per parity in Z')

# Consequence conditional on the T1D-Y & AA original-Weil analytic forms:
# positive selected E and L Schur blocks; AA positive conditioned Z sector.
# Sylvester inertia plus codimension -> negative index of original q on
# the fixed window smooth double-Mellin-pole-null core <=m<=2693201.
print('CONDITIONAL ANALYTIC COROLLARY: n_-(original q on D^0) <=',rank_upper)
print('NOT VERIFIED BY CODE: Fourier approximation, singular-value min-max,')
print('  original Weil source identities, prolate trace theorem, and full form-domain transfer')
print('OPEN: positivity on omitted low-frequency Z, complete RH Weil criterion, RH')
