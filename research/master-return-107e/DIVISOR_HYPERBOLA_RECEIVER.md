# RETURN107-E — Divisor hyperbola receiver with unconditional quadrature reserve

Status: new unconditional estimate for a component of the arithmetic
receiver, not a bound on the complete prime discrepancy. No numerical
search, native Lean proof, external referee acceptance or RH promotion.
The target is the literal tent of RETURN107-C/D, reopened before this
derivation. This uses only elementary convolution and quadrature.

## 1. Exact arithmetic split

Let X>=e^8, T(t)=(1-|t|)_+, and

    A(X)=sum_n Lambda(n)/sqrt(n) T(log(n/X)),
    c=int_R exp(t/2) T(t)dt =16 sinh(1/4)^2,
    c1=int_R t exp(t/2) T(t)dt.

The unconditional Dirichlet convolution identity Lambda=mu*log follows
from log=1*Lambda and mu*1=delta_1. All sums below are finite, because
dm lies in [X/e,eX]. Hence

    A(X)=sum_{d<=eX} mu(d)/sqrt(d)
                    sum_{m>=1} log(m)/sqrt(m) T(log(dm/X)).

Choose an integer cutoff D=floor X^(3/4). For d<=D put Y=X/d,
so Y>=X^(1/4)>=e². Define

    h_Y(t)=log(t)/sqrt(t) T(log(t/Y)), t>0,

extended by zero to t<=0. Its support is [Y/e,eY].
The continuous integral is exactly

    int_R h_Y(t)dt = sqrt(Y)[c log(Y)+c1].

The change of variable t=Y exp(u) gives this formula without a
Riemann approximation or a prime distribution assumption.

## 2. Unconditional lattice quadrature lemma

Claim: for Y>=e²,

    |sum_{m>=1} h_Y(m)-int h_Y(t)dt|
        <= 10 Y^(-3/2)(1+log Y).

Proof. The function h_Y is continuous, compactly supported, and
piecewise C². Its distributional second derivative is a finite
signed measure consisting of the ordinary h_Y'' and the jumps of
h_Y' at t=Y/e,Y,eY. On each open linear-log piece write
h_Y=f T(log(t/Y)) with f=t^-1/2 log t. Then

    f'=t^-3/2(1-(log t)/2),
    f''=t^-5/2((3/4)log t-2).

Since t>=Y/e>=e, |T|<=1, |dT/dt|<=1/t, and
|d²T/dt²|<=1/t² there,

    |h_Y''(t)|<=t^-5/2[(11/4)log t+4].

The jumps of h_Y' have magnitudes f(t)/t, 2 f(t)/t,
f(t)/t at the three breakpoints, respectively. Therefore the
total variation of the second derivative satisfies

    ||h_Y''||_TV
       <= e^(3/2) Y^(-3/2)
             [(35/6)log(eY)+8/3].

Here the integral contribution uses
int_(Y/e)^infty t^-5/2 dt=(2/3)(Y/e)^-3/2,
and the jump contribution uses 4(Y/e)^-3/2 log(eY).

On each unit cell [m,m+1], the exact trapezoid error is the
pairing of the second derivative measure with
(t-m)(m+1-t)/2, up to the conventional overall sign.
This kernel vanishes at both endpoints and has sup norm 1/8.
The identity follows by integration by parts twice for smooth
functions and extends to these piecewise C² continuous functions,
including derivative jumps, by the same integration by parts
or mollification. Summing over all integer cells gives

    |sum h_Y(m)-int h_Y| <= ||h_Y''||_TV/8.

The stated constant 10 is conservative: e^(3/2)<5,
(35/6)<6 and (35/6+8/3)=17/2<9, hence the preceding
coefficient is <= (5/8)(6 log Y+9)<10(1+log Y).
The sum over all integers equals the sum over positive integers
because Y/e>=e. This proves the lemma.

## 3. Receive the complete small-divisor contribution

Define the signed partial sums

    M0(D)=sum_{d<=D} mu(d)/d,
    M1(D)=sum_{d<=D} mu(d)log(d)/d.

Using the quadrature lemma, the small-divisor part A_small satisfies

    A_small(X)
      =sqrt X [(c log X+c1) M0(D)-c M1(D)] + R(X,D),

    |R(X,D)|<=10 X^(-3/2) sum_{d<=D} d(1+log(X/d))
      <=5 X^(-3/2) D(D+1)(1+log X).

For D=floor X^(3/4), D(D+1)<=2 X^(3/2), giving

    |R(X,D)|<=10(1+log X).

This is already within the RH-equivalent subpolynomial error
allowance: for every epsilon>0, 1+log X=O_epsilon(X^epsilon).
No cancellation of mu(d) was needed to achieve this reserve.

More generally D=floor X^theta with 0<theta<1 and X large
enough gives O(X^(2theta-3/2) log X). Thus the elementary
quadrature argument reaches theta=3/4 with only logarithmic
error; theta>3/4 is not justified at the same accuracy by
this particular total-variation bound. This is a limit of
the estimate, not a lower bound on the actual error.

## 4. The remaining signed hyperbola receiver

The large-divisor contribution is exactly

    H(X,D)=sum_{1<=m<=eX/D} log(m)/sqrt(m)
       sum_{d>D} mu(d)/sqrt(d) T(log(dm/X)).

The inner sum is finite, over d in [X/(em),eX/m] intersect
d>D; the m=1 term is zero. With D=floor X^(3/4), for X>=e^8
one has D>=X^(3/4)/2, so m<=2e X^(1/4).
Consequently

    A(X)-c sqrt X = S(X,D)+R(X,D),

    S(X,D)=sqrt X [(c log X+c1) M0(D)-c M1(D)-c]
               +H(X,D).

The preceding unconditional reserve proves that either of the
one-sided estimates

    S(X,floor X^(3/4)) <= O_epsilon(X^epsilon),
    S(X,floor X^(3/4)) >= -O_epsilon(X^epsilon)

for every epsilon>0, eventually, is sufficient for the matching
one-sided tent estimate of RETURN107-C, hence RH. Under RH the
tent discrepancy is bounded and the same decomposition gives
S=O(log X); thus the subpolynomial one-sided conditions are
also necessary. This is an exact reformulation with a genuinely
controlled remainder, not a proof of the signed bound.

Notice that the main signed combination and H must be received
together. Estimating each in absolute value, or replacing mu by
|mu|, destroys the relevant cancellation. For example from
|mu|<=1, summation over all d in the large-divisor window gives
an upper bound O(sqrt(X/m)) for the inner absolute sum (since
d is comparable to X/m); the outer bound becomes

    |H| <= O(sqrt X) sum_{m<=2e X^(1/4)} log(m)/m
          = O(sqrt X log² X).

This is at the wrong scale. Similarly |M0|<=1+log D,
|M1|<=O(log² D), yielding no improvement. These are explicit
bounds for the proposed elementary route, not an impossibility
theorem for all sieve or convolution methods. Positivity of
Lambda does not make the signed terms separately nonnegative.

## 5. Investigation decision

The quadrature error is paid unconditionally, at a cutoff
that leaves an outer sum of length O(X^(1/4)). The surviving
mathematical target is the cancellation between the truncated
Möbius Dirichlet moments and the signed large-divisor receiver.
No estimate sufficient to control that target has been
identified here. Importing RH-equivalent bounds on Möbius
sums as if they were unconditional is forbidden.

A subsequent investigation should either establish a new
signed cancellation inequality for S, with its dependence
on X explicit, or abandon this route if it only repackages
the original obstruction. More equivalent formulations
without an actual new inequality should not be counted as
proof progress. The fixed-aperture numerical evidence and
external mathematical review status are unchanged.

No claim of novelty for the elementary convolution or quadrature
method is made. No runner activated; RH remains unresolved.
