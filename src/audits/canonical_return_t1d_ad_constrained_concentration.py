#!/usr/bin/env python3
"""AD: exact rational checks only, conditional on AC/AB/AA/F/Y/Z analytic inputs."""
from fractions import Fraction as F
B=2**21; theta=F(1,2**38); N=1_380_000; a=F(241,120); k=14; C=3360; om=61000
assert F(314159,100000)*N/a-B>om
ratio=F(C,om)**(2*k)
U0=F(4*B,9*29*om)*ratio
U1=F(4*B,9)*ratio*(F(1+B*B,29*om)+F(2*B,28)+F(om,27))
assert U0/theta < F(1,2**80)
assert U1/theta < F(1,2**36)
t=F(1,2**40); d=F(1,2**18); eta=F(1,2**31)
# AC physical trace alpha(t,d)<10^-5
alpha=F(1,100000)
assert alpha>2*d and 4*(t*t/4+2*t*d)<(alpha-2*d)**2
# log2 < 7/10, exact exponential Taylor lower
term=s=F(1)
for n in range(1,7):
    term*=F(7,10*n); s+=term
assert s>2
assert 40*F(7,10)+F(1,2**40)<29
assert 82*F(7,10)+F(1,2**82)<58
assert 10*F(7,10)+F(1,2**10)<8
# log-energy and operator graph-norm of unprojected error
assert 94*t*t+F(10,2**40*100000**2)<F(1,10**20)
assert 94**2*t*t+F(94**2+2*94+2,3*2**40*100000**2)<F(1,10**18)
# orthogonal physical constraint carrier estimates
assert 2+2*37000**2+4*46**2<53000**2
assert B+53000<2**22
assert 2*9**2+2*193**2+4*11**2<275**2
assert F(1,10**10)+275*t<eta
# norm-one compressed low vector derivative <=2^41; alpha(1,2^41)<2^43
assert 4*(F(1,4)+2**42)<(2**43-2**42)**2
assert 65+58+F(2**86*59,3*2**82)<22**2
# operator norm of k(D) on sharp pole vectors, via log-squared tail
assert 4*(F(1,4)+2)<(5-2)**2 or 4*(F(1,4)+2)==(5-2)**2
assert 73**2+F(25*(73**2+2*73+2),3*2**10)<74**2
assert 2*74**2+2*37065**2+4*111**2<53000**2
M=1785+1073590
assert t*(2+t)<F(1,2**38)
assert 25*t*(2+t)<F(1,10**10)
assert 44*eta+eta**2<F(1,40_000_000)
assert 44*eta+eta**2+M*t*(2+t)<F(1,500_000)
assert F(1,10**9)+(53000+M)*t<F(1,950_000)
# Finite concentration split at theta/4 preserves the AA leakage condition.
assert F(1,2)+t<1
assert 2*N+1==2_760_001
assert F(1233,1750)>F(7,10)
print('CR-INT-T1D-AD EXACT-RATIONAL AUDIT PASS')
print('true constrained low: relative metric <2^-38; Weil energy <1/40000000')
print('relative second-Schur error <1/500000; mixed representer error <1/950000')
print('finite split rank <=2760001; high leakage <theta; no prolate angle required')
print('No spectral sign, mixed Gram entries or Lean certification claimed.')
