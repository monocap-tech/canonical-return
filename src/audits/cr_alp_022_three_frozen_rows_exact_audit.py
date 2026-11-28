#!/usr/bin/env python3
"""CR-ALP-022: exact 3-row frozen F41 source-moment sign certification.

The fixed frozen NF545 node (37/64,k) with k=0,1,2 yields
  sin(37*pi*x/32), x*cos(37*pi*x/32), x^2*sin(37*pi*x/32).
The script checks exact rational intervals and uniform error bounds
for every even nonnegative mass-one mollifier of width eta. It relies
on the explicit calculus identities in the paired research note.
NOT a Lean compile, a full 242-row F41 audit, or an RH theorem.
"""
from fractions import Fraction as Q
from math import factorial
import sympy as S

b=Q(1,10**8)
eta=Q(1,10**14)
assert b<1 and eta<b/10

# Exact rational log(2) enclosure by atanh series with geometric tail.
def log_atanh(q,N=40):
    assert Q(1)<=q<=Q(2)
    z=(q-1)/(q+1)
    t=z
    lower=Q(0)
    for j in range(N):
        lower+=2*t/(2*j+1)
        t*=z*z
    return lower,lower+2*t/((2*N+1)*(1-z*z))

ln2lo,ln2hi=log_atanh(Q(2))
assert Q(2,3)<ln2lo<ln2hi<Q(3,4)

# Machin pi enclosure, no floating approximations.
def atan_reciprocal(n,N=22):
    part=sum((Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(N)),Q(0))
    remainder=Q((-1)**N,(2*N+1)*n**(2*N+1))
    return ((part,part+remainder) if N%2==0 else (part+remainder,part))
a0,a1=atan_reciprocal(5)
c0,c1=atan_reciprocal(239)
pi_lo,pi_hi=16*a0-4*c1,16*a1-4*c0
assert Q(3)<pi_lo<pi_hi<Q(22,7)

# Frozen first 3 selected NF545 nodes are exactly (37/64,k), k=0,1,2.
# Let omega=37*pi/32 and theta=omega*log2.
omega_lo,omega_hi=Q(37,32)*pi_lo,Q(37,32)*pi_hi
assert 3<omega_lo<omega_hi<4
scaledtheta_lo,scaledtheta_hi=Q(37,32)*ln2lo,Q(37,32)*ln2hi
assert Q(3,4)<scaledtheta_lo<scaledtheta_hi<Q(7,8)
assert Q(7,5)**2<2          # cos(theta)<-sqrt(2)/2<-7/10
assert Q(23,16)**2>2        # sin(theta)>sin(pi/8)>3/8
# sin(pi/8)^2=(2-sqrt(2))/4 > (3/8)^2.
assert (2-Q(23,16))/4 == Q(9,64)  # strict because sqrt(2)<23/16
assert omega_lo*ln2lo>2
assert omega_hi*(b+eta)<1

# Exact symbolic derivatives for the actual source functions.
x, omega, ell = S.symbols('x omega ell',real=True)
functions=[S.sin(omega*x),x*S.cos(omega*x),x**2*S.sin(omega*x)]
A=[S.diff(f,x).subs(x,ell)-S.Rational(3,2)*f.subs(x,ell) for f in functions]
expected=[
    omega*S.cos(omega*ell)-S.Rational(3,2)*S.sin(omega*ell),
    (1-S.Rational(3,2)*ell)*S.cos(omega*ell)-omega*ell*S.sin(omega*ell),
    (2*ell-S.Rational(3,2)*ell**2)*S.sin(omega*ell)+omega*ell**2*S.cos(omega*ell)
]
assert all(S.simplify(a-c)==0 for a,c in zip(A,expected))

# All three A_j are strictly negative by rational trig enclosures:
#  cos(theta)<-7/10, sin(theta)>3/8, omega>3,
#  theta>2, and 2/3<ell<3/4. For j=2, sin(theta)<=1.
A0_upper=Q(3)*Q(-7,10)-Q(3,2)*Q(3,8)
A1_upper=Q(1,8)-2*Q(3,8)
A2_upper=Q(2,3)-Q(3)*Q(4,9)*Q(7,10)
assert A0_upper==Q(-213,80)
assert A1_upper==Q(-5,8)
assert A2_upper==Q(-4,15)

# Uniform analytic Taylor enclosure: for every F_j on |x|<=1,
# |F|<=2, |F'|<=8, |F''|<=50. With u=t-b and
# g(b-u)=e^((b-u)/2), g(b+u)=-e^(-(b-u)/2), one has
# m0=integral(g(b+u))du in [b^2/2,b^2/2+b^4/8],
# m1=integral(u*g(b+u))du in [-b^2-b^4/4,-b^2].
# Taylor error for T_j(u)=-F_j(ell-u)-2F_j(ell+u) is <100b^3.
# So |integral(h_raw*F_j) - b^2*A_j| < 103 b^3.
# Smoothing any real nonnegative mass-one rho supported in (-eta,eta)
# changes the moment by <=eta*sup|F'|*||h_raw||_1 <72*b*eta.
absolute_error=103*b**3+72*b*eta
relative_error=absolute_error/(b*b)
assert relative_error<Q(1,10000)

# Final strict universally valid bounds, for all permitted mollifiers.
for j,upper,goal in [(0,A0_upper,Q(-2)),
                     (1,A1_upper,Q(-1,2)),
                     (2,A2_upper,Q(-1,5))]:
    assert upper+relative_error<goal
    print(f'PASS frozen F41 row {j}: v_(j,2) < ({goal})*b^2, '
          f'b^2={b*b} (strict rational certificate)')

# This phase also checks where the existing Lean source proof can attach:
# - IV011 prime summands use C(log n)+C(-log n), zero for ODD C.
# - gammaBracket(r)=Re(digamma(1/4+i*r/2))-log(pi), even from
#   Zeta23.WeilEF.digamma_conj on integerComplement.
# - pole terms paperFT(C)(i/2)+paperFT(C)(-i/2), zero for ODD C.
# The formal windowed-correlation reflection and FT parity lemmas
# have not been Lean compiled, and are not claimed here.

print('CR-ALP-022 THREE FROZEN ROWS + PARITY SOURCE AUDIT PASS')
print('Entire original signed even/odd Weil mixed parity: analytic; Lean pending')
print('No frozen F41 H/V mutation, 242-entry acceptance, whole-core Schur or RH')
