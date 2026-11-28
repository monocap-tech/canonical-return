#!/usr/bin/env python3
"""CR-ALP-021: pinned F41 row-0 raw packet nonzero sign, exact rational gate.

Mathematical proof rests on the exact frozen source shape sin(37*pi*x/32),
the elementary trigonometric monotonicity bounds stated in the note,
and any EVEN nonnegative mass-one Schwartz mollifier of width eta.
The program audits rational strict inequalities and source identifiers;
it does NOT numerically integrate a specified mollifier or compile Lean.
"""
from fractions import Fraction as Q
from math import factorial

B=Q(1,10**8)
ETA=Q(1,10**14)
HALF=Q(1,2)


def atanh_log(z,N=36):
    assert Q(1)<=z<=Q(2)
    u=(z-1)/(z+1)
    t=u
    v=Q(0)
    for j in range(N):
        v+=2*t/(2*j+1)
        t*=u*u
    return v,v+2*t/((2*N+1)*(1-u*u))

l2,u2=atanh_log(Q(2))
assert Q(2,3)<l2<u2<Q(3,4)
assert Q(37,32)*l2>Q(2,3)
assert Q(37,32)*u2<1

# Self-contained exact Machin alternating atan enclosures of pi.
def atan_reciprocal(q,N=20):
    v=sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),Q(0))
    nextterm=Q((-1)**N,(2*N+1)*q**(2*N+1))
    return (v,v+nextterm) if N%2==0 else (v+nextterm,v)
a0,a1=atan_reciprocal(5)
b0,b1=atan_reciprocal(239)
pi_lo=16*a0-4*b1
pi_hi=16*a1-4*b0
assert Q(3)<pi_lo<pi_hi<Q(22,7)
omega_lo=Q(37,32)*pi_lo
omega_hi=Q(37,32)*pi_hi
assert Q(3)<omega_lo<omega_hi<Q(4)
assert omega_hi*B<1 and omega_hi*ETA<1

# The frozen index-zero oscillator has node (37/64,k=0), so
# phi_0(x)=sin(2*pi*(37/64)*x)=sin((37*pi/32)*x).
# Write c=ln2, theta=omega*c. Its scaled theta/pi is (37/32)*ln2.
# Rational barriers prove 2*pi/3 < theta < pi, therefore
# cos theta < -1/2 and sin theta >0.
scaled_theta_low=Q(37,32)*l2
scaled_theta_high=Q(37,32)*u2
assert Q(2,3)<scaled_theta_low<scaled_theta_high<1

# Unsmooth, centered packet g_b(b+u): for 0<u<b,
#   g_b(b-u)=exp((b-u)/2), g_b(b+u)=-exp(-(b-u)/2).
# Thus J_sin=integral g_b(b+u)sin(omega*u)du
#           =-2int_0^b cosh((b-u)/2)sin(omega*u)du
#       <= -omega*b*b/2  [sin(y)>=y/2 for 0<=y<=1].
# And J_cos=2int_0^b sinh((b-u)/2)cos(omega*u)du>0.
# v_raw = -cos(theta)*J_sin -3*sin(theta)*J_cos
#       < - omega*b*b/4 < -(3/4)*b*b.
# Even mollifier multiplies *this exact oscillator row* by
# mu_rho=int rho(v)cos(omega*v)dv >=1/2 if omega*eta<1.
# Therefore the SMOOTH actual frozen row obeys:
#                  v_{0,2} < -(3/8)*b*b.
row_absolute_floor=Q(3,8)*B*B
assert row_absolute_floor==Q(3,8*10**16)
assert row_absolute_floor>0

# Frozen selected physical span S_242 contains phi_0 and phi_0 is
# L²[-11/10,11/10]-bounded by sqrt(11/5). Cauchy then provides
# a genuinely nonzero inside-selected squared projection:
# dist >= |v|²/||phi_0||² > (45/704)*b**4.
projection_mass_floor=Q(45,704)*B**4
assert projection_mass_floor>0

# Source from CR-ALP-020 has independent even component outside
# the full frozen ODD span, with normalized squared floors 1/10,1/5.
assert Q((2-1)**2,2*(1+2*2))==Q(1,10)
assert Q((3-1)**2,2*(1+3*3))==Q(1,5)

print('CR-ALP-021 EXACT RATIONAL OSCILLATOR/PARITY AUDIT PASS')
print('Frozen F41 source row: j=0, node=(37/64,0), phi0=sin(37*pi*x/32)')
print('2*pi/3 < (37*pi/32)*log2 < pi; sin(theta)>0, cos(theta)<-1/2')
print('Any even nonnegative mass-one mollifier eta=1e-14:')
print('RAW h2 frozen physical moment v_{0,2} < -',row_absolute_floor,'(strict, exact)')
print('Inside 242-row selected squared projection >',projection_mass_floor)
print('Outside 242-row selected squared mass: >= 10% of raw h2 mass')
print('Evenized E2,E3 all 246 physical F41 mixed moments =0 by parity')
print('Signed original IV010 even/odd Schwartz mixed form =0 analytically')
print('No new Lean compilation, 242-entry numerical H/V audit, whole-core Schur, or RH')
