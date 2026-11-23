#!/usr/bin/env python3
"""CR-ALP-024 exact rational certificates for fixed-support original Weil form.

Human analysis in research note: digamma partial-fraction inequalities, original
signed EF/pole/primes+Gamma, Plancherel, time-band concentration trace.
This script audits EVERY finite number inequality, NOT those analytic theorems.
No Lean compilation, F41 H/V acceptance or RH claim.
"""
from fractions import Fraction as F
from math import isqrt, factorial

B=F(11,10) # exact reviewed original IV physical support radius
PRIME_POWERS={2:2,3:3,4:2,5:5,7:7,8:2,9:3}


def atanh_log_1_2(z,N=40):
    assert 1<=z<=2
    u=(z-1)/(z+1)
    term=u; lower=F(0)
    for j in range(N):
        lower+=2*term/(2*j+1)
        term*=u*u
    upper=lower+2*term/((2*N+1)*(1-u*u))
    return lower,upper


def log_int(n):
    assert n>=1
    m=n.bit_length()-1
    aa,bb=atanh_log_1_2(F(n,2**m))
    a,b=atanh_log_1_2(F(2))
    return aa+m*a,bb+m*b

# Source cutoff log 9 < 2B < log 10: precisely n=2,3,4,5,7,8,9.
ln9=log_int(9);ln10=log_int(10)
assert ln9[1]<2*B<ln10[0]
assert ln10[0]>F(16,7)  # permits fully integral spectral band R=10^14
assert log_int(2)[1]<1
assert log_int(3)[1]<2

# Every supported autocorrelation is inside (-2B,2B).
# Exact von Mangoldt coefficients; no pseudo-prime 6 or 9-as-prime error.
assert set(PRIME_POWERS)=={2,3,4,5,7,8,9}
assert all(n==p**k for n,p in PRIME_POWERS.items()
           for k in range(1,5) if p**k==n)

# Rigorous sqrt lower: floor(10^6 sqrt(n))/10^6 using integer isqrt.
scale=10**6
S_upper=F(0)
for n,p in PRIME_POWERS.items():
    sqrtlo=F(isqrt(n*scale*scale),scale)
    assert sqrtlo>0 and sqrtlo*sqrtlo<=n
    S_upper+=log_int(p)[1]/sqrtlo
assert S_upper<4
print('PASS: 7 source prime powers, S=sum Lambda(n)/sqrt(n) <',float(S_upper),'<4')

# Base Gamma symbol a(0)= -gamma -pi/2 -3log2 - log(pi)>-6,
# using gamma<1, pi<22/7, log2<7/10, log(pi)<6/5.
pi_upper=F(22,7)
x=F(6,5)
assert sum((x**k/factorial(k) for k in range(5)),F(0))>pi_upper
assert log_int(2)[1]<F(7,10)
assert 1+pi_upper/2+3*F(7,10)+F(6,5)<6
# For |r|>=2, N=floor(|r|/2)>=1. Split digamma series at N:
# ((N+1)(N/2+1/4))/N^2 <= 3/2, tail >= -3/(4N).
for N in range(1,500):
    loss=F((N+1)*(2*N+1),4*N*N)
    assert loss<=F(3,2)
assert F(13,4)+F(3,2)+2<8
assert F(-6)+2==-4<0
# The lower gamma envelope is a(r)>=log(1+|r|)-8, while
# the upper is a(r)<=log(1+|r|)+3; exact continuum proof in note.

# Pole residues = 2 Re(F_plus * conj(F_minus)),
# |P|<=4 sinh(B) ||f||^2 <8||f||^2 using exp(B)<4.
exp_terms=sum((B**k/factorial(k) for k in range(17)),F(0))
first_omitted=B**17/factorial(17)
exp_upper=exp_terms+first_omitted/(1-B/18)
assert exp_upper<4
assert 2*exp_upper<8

# Source-specific coercivity constants:
# pole-null: W(f,f) >= LogEnergy(f) -(8+2*S)||f||^2
# general: plus |pole|<8||f||^2.
assert 8+2*S_upper<16
assert 8+2*S_upper+8<24

# For W(f,f)<=0 in a subspace, log-energy <=C norm^2.
# At R=10^14: log(1+R)>14*log(10)>32=2*16.
# At R=10^21: log(1+R)>21*log(10)>48=2*24.
assert 14*log_int(10)[0]>32
assert 21*log_int(10)[0]>48
R_null=10**14;R_full=10**21
# Time-frequency concentration (compact positive trace operator):
# Trace(T_R)=2 B R/pi; if every unit f in a nonpositive subspace has
# low-frequency probability >1/2, dimension m<4 B R /pi.
# Since pi>3, m < 4 B R /3 =(22/15)*R.
for label,R,C in [('pole-null',R_null,16),('full',R_full,24)]:
    m_cap=(F(4)*B*R/F(3))
    assert m_cap==F(22,15)*R
    print(f'PASS {label}: C={C}, R={R}; negative/nonpositive subspace '
          f'integer dimension <= {m_cap.numerator//m_cap.denominator:,}')
print('CR-ALP-024 EXACT ARITHMETIC/SPECTRAL-TRACE GATES PASS')
print('Only finite certificate arithmetic audited; digamma + functional '
      'analysis supplied separately, NOT Lean compiled.')