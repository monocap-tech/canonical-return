# RETURN107-D — Exact Riesz primitive, growth obstruction and prime-power separation

Status: internal analytic derivation; no unconditional RH-strength growth bound
and no theorem promotion. Uses RETURN107-C's literal Q and tent receiver.
The purpose is to determine precisely what an arithmetic control argument
must accomplish, without relying on a formal zero sum or a contour shift.

## 1. Absolutely convergent centered partial fractions

Let z=s-1/2, a_rho=rho-1/2 for each nontrivial zero counted with
multiplicity, and beta_m=2m+1/2 for m>=1. Define

    F(z)=-zeta'(1/2+z)/zeta(1/2+z).

F is analytic at zero: zeta(1/2)!=0 follows from its alternating eta
series. Use the classical Hadamard product for the order-one completed
xi function and the digamma partial fractions. These are standard
unconditional zeta identities. In particular the difference of the
xi logarithmic derivative at s=1/2+z and at s=1/2 equals

    sum_rho [1/(z-a_rho)+1/a_rho].

The difference converges locally normally away from the poles since
sum |a_rho|^-2<infty by the classical N(T)=O(T log T) zero-counting bound.
Subtracting the affine Taylor polynomial of the logarithmic derivative
at z=0 is therefore legitimate termwise. For any nonzero pole location a,

    1/(z-a) + 1/a + z/a² = z²/[a²(z-a)].

The identity
F(z)=-xi'/xi(1/2+z)+1/(z+1/2)+1/(z-1/2)
       -(log pi)/2+psi(1/4+z/2)/2
has a cancellation of the pole at z=-1/2 between the first rational
term and the m=0 digamma term. The remaining rational pole is
z=1/2 with residue +1; all nontrivial and trivial zero poles have
residue minus their multiplicity. Thus

    F(z)/z² = F(0)/z² + F'(0)/z + 4/(z-1/2)
               -sum_rho 1/[a_rho²(z-a_rho)]
               -sum_{m>=1} 1/[beta_m²(z+beta_m)].

The last two sums are locally normally convergent away from their
poles. No enumeration of zero positions is necessary. At the center
xi'/xi(1/2)=0 by its symmetry, giving F(0)=g(0)/2, where
g(0)=psi(1/4)-log pi is the physical Gamma symbol constant.

## 2. Exact primitive via Laplace uniqueness

For Re z>1/2, Euler product and absolute convergence give

    integral_0^infty exp(-zL)
         sum_{log n<=L} Lambda(n)/sqrt(n)(L-log n) dL
      = F(z)/z².

Indeed the individual ramp transform is n^-z/z². The exponential
main term 4 exp(L/2) has transform 4/(z-1/2). Define

    Q(L)=sum_{log n<=L} Lambda(n)/sqrt(n)(L-log n)-4 exp(L/2),
    q(L)=Q(L)-F(0)L-F'(0).

The partial fraction identity proves

    q(L)=-sum_rho exp(a_rho L)/a_rho²
         -sum_{m>=1} exp(-beta_m L)/beta_m²,  L>=0.

To justify inversion without a contour shift, the series on the
right is absolutely uniformly convergent on every compact L interval
by 0<Re rho<1 and sum |a_rho|^-2<infty. It is of exponential order
1/2 on L>=0. Its Laplace transform can be taken termwise on
Re z>1/2, by summable absolute majorants. Q is continuous and its
Laplace transform there is already justified by the absolutely
convergent ramp integrals. Laplace uniqueness for continuous
functions of exponential order (or Fourier uniqueness on an
exponentially damped L1 difference) yields the identity everywhere.
The exact affine constant is not arbitrary: at L=0 the original
Q(0)=-4 fixes the identity F'(0)=-4+sum_rho a_rho^-2+sum beta_m^-2.
No RH assumption is used.

The trivial-zero term is bounded in absolute value by
exp(-5L/2) sum beta_m^-2 for L>=0. The nontrivial series is real
because zeros occur in conjugate pairs. Its absolute convergence
does not imply boundedness uniformly on the half-line: the
factor exp(Re(a_rho)L) cannot be dropped.

## 3. Exact growth frontier of the centered primitive

Set Theta=sup{Re rho: rho a nontrivial zero}, and sigma=Theta-1/2.
By the functional equation sigma>=0 and sigma<=1/2. The absolute
series gives the unconditional bound

    |q(L)| <= (sum |a_rho|^-2 + sum beta_m^-2) exp(sigma L).

Under RH sigma=0, so q is bounded. Conversely if there is a zero
with Re a_rho>delta>=0, a one-sided eventual estimate

    q(L)<=M exp(delta L), or q(L)>=-M exp(delta L)

is impossible. Here M>=0 can be chosen without loss. To prove it,
choose L0>0 beyond the eventual bound, and take the nonnegative
function M exp(delta L) minus or plus q(L). The Laplace transform
of q on [L0,infty) is the transform on [0,infty) minus an entire
compact-interval integral. The partial fractions show a nonzero
pole at each such a_rho, of residue -m_rho/a_rho², with the
positive real pole 1/2 already subtracted. Hence the abscissa of
convergence alpha of the nonnegative transform must be at least
Re a_rho>delta. But the transform is holomorphic near real alpha:
zeta has no positive real centered zero and the pole of the
added exponential is at delta. The Landau lemma proved in
RETURN107-B forbids this. Therefore, if sigma>0, for every
0<=delta<sigma,

    limsup_{L->infty} exp(-delta L) q(L)=+infty,
    liminf_{L->infty} exp(-delta L) q(L)=-infty.

This is a conditional obstruction, not a detection of an actual
off-critical zero. It proves that either one-sided eventual
boundedness of q, or one-sided O(exp(epsilon L)) for every
epsilon>0, is equivalent to RH. Even either one-sided O(1) bound
for the original Q after removing only F(0)L is equivalent.
The affine term is immaterial at any strictly positive exponent.

The exact tent curvature has the same zero obstruction: applying
the step-b second difference divided by b² to the series is
legitimate for L>=b, and its multiplier is
(exp(b a)-2+exp(-b a))/b². At a nontrivial pole it is nonzero
when Re a!=0; the centered affine part cancels. This confirms
that the cumulative and tent receivers have the same exponential
frontier, not that their numerical norms are equal.

## 4. Isolate the actual primes in the elementary tent

Let T(t)=(1-|t|)_+, X>=e, L=log X, and

    A(X)=sum_n Lambda(n)/sqrt(n) T(log(n/X)),
    c_T=16 sinh(1/4)^2.

Write A=A_1+A_2+A_ge3 by n=p^m and m=1,2,>=3.
All terms nonnegative. For m>=3 contributing n in [X/e,eX],
the exponent is <=log(eX)/log 2 and p<=(eX)^(1/3).
Each weight is <=sqrt(e/X) log(eX). Counting all integers
instead of primes gives the unconditional crude bound

    0<=A_ge3(X)=O(X^(-1/6) log² X)=o(1).

For squares put theta(t)=sum_{p<=t} log p and write

    A_2(X)=integral_{sqrt(X/e)}^{sqrt(eX)}
             T(2log t-L)/t d theta(t).

The classical unconditional PNT theta(t)=t+o(t) gives
A_2(X)=1/2+o(1). To see the error explicitly, the test
w(t)=T(2log t-L)/t is continuous, piecewise C1, vanishes
at the integration endpoints, and |w'(t)|<=3/t² away
from its finite breakpoints. Uniformly on this fixed ratio
interval, |theta(t)-t|=o(sqrt X). Stieltjes integration by
parts therefore bounds the difference from the dt integral
by o(sqrt X)*O(X^-1/2)=o(1). The dt integral is exactly

    integral T(2log t-L) dt/t = (1/2) integral T(u)du = 1/2.

This uses only the ordinary PNT, not RH or a square-root
error estimate. Hence

    A(X)=A_1(X)+1/2+o(1),

where A_1 uses primes only. For example the upper one-sided
RH-equivalent condition becomes

    sum_{X/e<=p<=eX} log(p)/sqrt(p) T(log(p/X))
         <= c_T sqrt(X)-1/2+O(1),

and likewise for the lower side or O_epsilon(X^epsilon).
The displayed -1/2 is a genuine asymptotic contribution,
although a bounded-error criterion is insensitive to any
fixed shift. No estimate of the required prime-only error
has been established here.

Even without invoking PNT for the limiting square offset,
Chebyshev's elementary theta(t)=O(t) and partial summation
give A_2=O(1), sufficient to discard all higher powers for
a bounded-error equivalence. Thus the RH-level obligation
does not reside in prime-power corrections.

## 5. Investigation decision and boundary

The centered Riesz expansion is exact, and the growth obstruction
can be stated as a conditional two-sided exponential oscillation
result. The ordinary PNT pays the entire square correction, and
higher powers decay without PNT. The residual arithmetic problem
is a one-sided estimate for the smoothed distribution of primes
themselves. No use of the zero expansion as an unconditional
bound is allowed: doing so would require exclusion of exactly
the off-critical exponent that is in question.

Next investigate whether a proposed arithmetic estimate can
control this prime-only discrepancy without circular use of an
RH-strength zero-free assertion. A new equivalent formulation
alone is not progress toward that estimate. No novelty relative
to classical explicit-formula and oscillation results is claimed;
the specific goal here is compatibility with RETURN's literal
receiving object. No new runner, local proof promotion, or RH
resolution.
