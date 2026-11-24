#!/usr/bin/env python3
"""CR-INT-T1D-AC exact-fraction arithmetic audit (not a proof of analytic imports)."""
from fractions import Fraction as F

B, theta, a = 2**21, F(1, 2**38), F(241,120)
N, k, C, omega_floor = 1_380_000, 14, 3360, 61000
pi_low = F(314159,100000)  # imported T1D-AB Machin enclosure
assert pi_low*N/a-B > omega_floor

# T1D-AB coefficient tail and H1-weighted Fourier-tail refinement.
ratio = F(C,omega_floor)**(2*k)
e0sq = F(4*B, 9*(2*k+1)*omega_floor)*ratio
e1sq = F(4*B,9)*ratio*(
    F(1+B*B, (2*k+1)*omega_floor)
    +F(2*B,2*k)+F(omega_floor,2*k-1))
assert e0sq/theta < F(1,2**80)
assert e1sq/theta < F(1,2**36)
t, d = F(1,2**40), F(1,2**18)

# Zero-extension boundary Fourier tail:
# alpha=2sqrt(t*t/4+2*t*d)+2*d < 1/100000.
alpha_cap = F(1,100_000)
inner = t*t/4 + 2*t*d
assert alpha_cap-2*d > 0
assert 4*inner < (alpha_cap-2*d)**2

# Positive rational Taylor sum for exp(0.7)>2, hence log(2)<0.7.
x=F(7,10); s=F(0); term=F(1)
for n in range(7):
    if n: term=term*x/n
    s+=term
assert s>2

# log(1+2^40)<29, log(1+2^-40)<2^-40.
Rz=2**40
Eerr_sq_upper=94*t*t+F(10, Rz*100_000**2)
assert Eerr_sq_upper < F(1,10**20)   # energy-norm error <10^-10

# For e in the excluded prolate eigenspace:
# ||e'|| <= B theta^-1/2 ||e|| = 2^40||e||;
# alpha(1,2^40)<2^42 and log(1+2^80)<57.
Rbig=2**80
Ebase_sq_upper=65+57+F(2**84*58,3*Rbig)
assert Ebase_sq_upper < 22**2

# Normalized physical carrier log-energy norms:
# two pole carriers <9, two selected <193, four central/satellite <11.
assert F(65+8)+F(25*9, 3*2**10) < 9**2
assert 65+37000 < 193**2
assert 65+46 < 11**2
projection_carrier_energy_sq=2*9**2+2*193**2+4*11**2
assert projection_carrier_energy_sq < 275**2
projected_energy_err=F(1,10**10)+275*t
assert projected_energy_err < F(1,2**31)
assert 22+275 == 297

# Full original Gamma+signed-prime form (on the pole-null projection).
q_diff_upper=F(594,2**31)+F(1,2**62)
assert q_diff_upper < F(1,3_000_000)
metric_diff_upper=t*(2+t)
assert metric_diff_upper < F(1,2**38)
assert 25*metric_diff_upper < F(1,10**10)

# Imported global R,T second-Schur norms from T1D-F/Z.
R_schur_op,T_schur_op=1785,1_073_590
assert 2+t < F(201,100)
second_schur_form_diff=F(1,3_000_000) + (
    R_schur_op+T_schur_op)*F(201,100)*t
assert second_schur_form_diff < F(1,400_000)

print("CR-INT-T1D-AC EXACT-RATIONAL AUDIT PASS")
print("finite frame size <=",2*N+1)
print("L2 error squared <2^-80; derivative error squared <2^-36")
print("physical metric squared-norm error <2^-38 ||e||^2")
print("signed prime error <10^-10 ||e||^2")
print("complete Weil/Gamma error <1/3000000 ||e||^2")
print("complete second-Schur error <1/400000 ||e||^2")
print("Unproved: prolate frame angle, low/high mixed Schur, spectral sign.")
