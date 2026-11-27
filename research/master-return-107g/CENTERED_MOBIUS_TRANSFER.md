# RETURN107-G — Centered Möbius receiver and weighted transfer estimate

Status: analytic component estimate and explicit route assessment.
No new runner, native Lean proof, external referee acceptance, or RH
promotion. Reopened RETURN107-F and use exactly its fixed smooth
dyadic autocorrelation C, supported in [-2,2].

## 1. Center the entire divisor receiver, not just a cutoff

Let L=log X, c=int exp(s/2) C(s)ds and
c1=int s exp(s/2) C(s)ds as in RETURN107-F. Define for all real v

    r(v)=exp(-v/2) sum_{m>=1} log(m)/sqrt(m) C(log m-v)
            -cv-c1.

The sum is locally finite, hence r is smooth. For v<=-2 the
lattice sum is zero (the m=1 term is zero), so r(v)=-cv-c1.
For v>=20, RETURN107-F gives

    |r(v)|<=2 exp(5)(1+v) exp(-v²/20).

The normalized receiver has the exact representation

    A_C(X)/sqrt X
      =sum_{d<=e² X} mu(d)/d
           [c(L-log d)+c1+r(L-log d)].

The bracket vanishes for d>e² X, so it can be extended to all d
as a finite sum. To separate its terms requires convergence,
which we establish using an explicitly external classical estimate:

    M(t)=sum_{n<=t} mu(n)=O(t exp(-k sqrt(log t))) for some k>0.

This is the classical unconditional Mertens-function bound from
the zeta zero-free region, not RH. See the statement and discussion
in "On some effective results involving zeros of the Riemann zeta
function", Bulletin of the Australian Mathematical Society,
https://www.cambridge.org/core/journals/bulletin-of-the-australian-mathematical-society/article/on-some-effective-results-involving-zeros-of-the-riemann-zeta-function/8D885EAF542CEA62250CDF37719DC7B0 .
We use the weaker classical square-root logarithmic estimate,
not a claim of a new best known Mertens bound.

Partial summation gives convergence of sum mu(d)/d and
sum mu(d) log(d)/d. In fact the integrated tails are dominated
by integrals of (1+log t) exp(-k sqrt(log t)) dt/t, finite.
The same majorant, uniform for real s>=1 near 1, allows passage
to the limit in 1/zeta(s) and its derivative. Since zeta has
a simple pole of residue one at 1, the exact identities are

    sum mu(d)/d=0,
    sum mu(d)log(d)/d=-1.

For clarity, the derivative is -sum mu(d)log(d)/d, and its
limit at s=1 is +1. The sums are in natural integer order,
not absolutely convergent rearrangements. The r-series
converges by the same tail argument, because r(L-log d)
is affine in log d for sufficiently large d.

Thus the complete signed combination, including the previous
truncated Dirichlet moments, can be centered exactly:

    E_C(X):=A_C(X)-c sqrt X
          =sqrt X sum_{d>=1} mu(d)/d r(log(X/d)).

No main term or large-divisor correction remains outside this
single response weight. It is a conditional-convergence identity
with the classical unconditional bound as an explicit premise.

## 2. Derivative reserve and exact integral transport

Set K(v)=r(v)+r'(v). For v<=-2,

    K(v)=-cv-c1-c.

We need quantitative decay on the opposite side, including the
derivative rather than just r itself. Differentiating the
locally finite normalized lattice sum shows r' is the same
normalized quadrature error with C replaced by
C_d=-C'/1-C/2 (that is, -C'-C/2).
Its continuous constants are

    int exp(s/2) C_d(s)ds=0,
    int s exp(s/2) C_d(s)ds=c,

by integration by parts (C is smooth compactly supported).
Therefore the derivative of the continuous term is exactly c,
as required. For any fixed integer N>=2, applying the scaling,
chain-rule and Poisson argument of RETURN107-F to C_d gives

    |r'(v)|<=B_N (1+v) exp(-Nv), v>=0.

Here B_N is finite because every derivative through N of C_d
is integrable. More explicitly use the polynomial
prod_{j=0}^{N-1}(d/ds-j-1/2) and the same L1 estimate,
with max_{k<=N} ||C_d^(k)||_1 in place of the earlier bound
for C. Leibniz adds at most N times the same maximum.
The identical estimate holds for r by the same argument.
Thus for some finite B (choose N=2),

    |K(v)|<=B(1+v) exp(-2v), v>=0,
    |K(v)|<=B(1+|v|), v<=0.

The second inequality follows from the affine formula below -2
and continuity on [-2,0]. These bounds do not require a
numerical evaluation of C or its derivatives.

Abel summation of the centered series now gives exactly

    E_C(X)=sqrt X int_1^infty M(t)/t² K(log(X/t)) dt.

Indeed the derivative of t^-1 r(log(X/t)) is
-t^-2 [r(log(X/t))+r'(log(X/t))]. The lower endpoint uses
M(1^-)=0, and the upper boundary term vanishes by the
classical bound since r is only logarithmically growing at
infinity in t. The integral is absolutely convergent with
that bound. The transfer is valid for every X>=1.
Equivalently, changing variables v=log(X/t),

    E_C(X)=sqrt X int_{-infty}^{L} [M(X exp(-v))/(X exp(-v))]
                                                 K(v) dv.

The sign of K is not assumed, and the integral is not asserted
to be a positivity-preserving operator. In particular the
negative-v tail extends indefinitely; replacing it by a compact
window would omit the very signed Dirichlet-moment contribution
being centered.

## 3. Quantitative unconditional transfer

Write M(t)/t <= in absolute value B0 exp(-k sqrt(log t))
for t>=1, enlarging B0 for bounded t. Split the v-integral at 0.
For 0<=v<=L, sqrt L-sqrt(L-v)<=sqrt v, so

    exp(-k sqrt(L-v)) <= exp(-k sqrt L) exp(k sqrt v).

The positive-v part is bounded by a finite constant times
exp(-k sqrt L), using int_0^infty (1+v) exp(-2v+k sqrt v)dv.
For v=-w<=0, the negative-v part is bounded by

    B B0 int_0^infty (1+w) exp(-k sqrt(L+w)) dw.

Substitute y=sqrt(L+w)-sqrt L. Then
w=2 sqrt L y+y², dw=2(sqrt L+y)dy. The resulting polynomial
in sqrt L and y times exp(-ky) integrates to
O_k((1+L) exp(-k sqrt L)). All constants depend only on C
and the fixed unconditional Mertens bound. Hence

    |E_C(X)|=O(sqrt X (1+log X) exp(-k sqrt(log X))).

This is a genuine unconditional bound for the full signed
receiver but follows from an already-known arithmetic theorem.
It is NOT an improvement on the classical prime number theorem
error estimates, nor is it subpolynomial in X. In particular
its exponent of X tends to 1/2, not zero.

## 4. Weighted power transport without logarithmic loss

For 0<=theta<1 define the finite positive constant

    J_theta=int_R exp((1-theta)v) |K(v)| dv.

It is finite: on v>=0 the bound is O((1+v)exp(-(1+theta)v));
on v<=0 it is O((1+|v|)exp((1-theta)v)).
If an independent arithmetic estimate gives
|M(t)|<=B_theta t^theta for every t>=1, the exact transfer yields

    |E_C(X)|<= B_theta J_theta X^(theta-1/2).

This avoids the additional logarithmic factors that arise
from bounding the truncated moments and large-divisor term
separately. It gives a precise weighted-norm transfer, not
a new estimate for M(t). In particular the hypothesis
M(t)=O_epsilon(t^(1/2+epsilon)) for all epsilon>0 is itself
RH-strength: the series for 1/zeta(s), continued by Abel
summation, would be holomorphic in Re s>1/2, excluding zeros
there by analytic uniqueness. Do not assume it as a
routine arithmetic input. Conversely the known RH implication
for M is not needed for this receiver's sufficient direction.

## 5. Assessment of the short-outer-sum route

RETURN107-F paid the small-divisor lattice error with o(1),
and shortened the remaining outer index range to X^o(1).
This stage pays the centering of the truncated signed moments
and the long inner windows exactly by one Abel integral.
It provides an unconditional full-receiver estimate and a
weighted transfer estimate with no logarithmic loss, but
neither the short outer range nor the absolute kernel estimate
supplies the missing RH-strength cancellation. The negative-v
tail of K is affine, not rapidly decaying before the Mertens
weight is applied. No uniform local-to-global contraction has
been established.

Decision: stop treating the divisor-quadrature route as a
closing mechanism unless a genuinely new signed cancellation
estimate is supplied. Further routine cutoff optimization or
reformulation of the centered receiver is not evidence of
progress toward RH. A subsequent stage should compare this
global obstruction with another existing transport route,
rather than silently assuming the Mertens estimate required
to close it. The fixed-aperture numerical acceptance package
is unaffected. No novelty claim for the classical arithmetic
input or Abel summation. RH remains unresolved.
