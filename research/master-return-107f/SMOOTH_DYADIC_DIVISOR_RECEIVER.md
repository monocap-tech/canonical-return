# RETURN107-F — Smooth dyadic divisor receiver with subpolynomial outer range

Status: unconditional analytic quadrature refinement of RETURN107-E,
not a bound on the signed arithmetic discrepancy. No new runner,
native Lean certificate or RH promotion. The exact smooth receiver
of RETURN107-B is reopened; this changes the kernel from the tent,
so no identity of their arithmetic discrepancies is asserted.

## 1. Receiver and derivative bounds

Use b=1, a_j=2^-j, v_j=1_{[-a_j,a_j]}/(2a_j), and the smooth
probability density u_*=v_1*v_2*... of RETURN107-B. Set
C=u_* * u_*, using evenness. It has support [-2,2], mass one,
and bilateral Laplace transform U(z)^2, nonzero for Re z!=0.
Its one-sided subpolynomial prime-discrepancy criterion is
therefore already attached to RH in RETURN107-B.

For every integer k>=0,

    ||C^(k)||_1 <= 2^(k(k+1)/2).

Proof: distribute one derivative to each of the first k distinct
uniform factors in one copy of u_*. The distributional derivative
of v_j is the difference of endpoint point masses divided by
2a_j, of total variation 1/a_j=2^j. The remaining infinite
convolution is a probability measure by the convergent random
series construction, as is the other copy of u_*. Convolution
of finite signed measures has total variation at most the product
of their total variations. The derivative identity follows by
distributional convolution of compactly supported measures.
Since C is smooth, its derivative measure has the indicated
L1 norm. For k=0 use mass one.

For Y>=e^20 define h_Y(t)=t^(-1/2) log(t) C(log(t/Y))
for t>0, extended by zero on t<=0. It is smooth compactly
supported in [Y exp(-2),Y exp(2)]. Write

    h_Y(t)=Y^(-1/2) g_Y(t/Y),
    g_Y(v)=v^(-1/2)(log Y+log v) C(log v).

Let s=log v and P_Y(s)=(log Y+s) C(s). The exact derivative
formula for N>=1 is

    g_Y^(N)(v)=exp(-(N+1/2)s)
               product_{j=0}^{N-1}(d/ds-j-1/2) P_Y(s).

It follows by induction using d/dv=exp(-s)d/ds.
Changing variables in the L1 norm gives a prefactor
exp((1/2-N)s), bounded by exp(2N+1) on [-2,2].
The sum of absolute polynomial coefficients of the indicated
product is product_{j=0}^{N-1}(j+3/2) <=(N+1)^N.
For 0<=k<=N, Leibniz's rule and the derivative bound give

    ||P_Y^(k)||_1 <= (log Y+2+N) 2^(N(N+1)/2).

Consequently

    ||g_Y^(N)||_1
      <= exp(2N+1)(N+1)^N 2^(N(N+1)/2)(log Y+2+N)
      <= exp(5N²)(log Y+2+N).

The last inequality uses 2N+1<=3N², log(N+1)<=N,
and log(2) N(N+1)/2<=N² for N>=1. The original scale is

    ||h_Y^(N)||_1 <= Y^(1/2-N) exp(5N²)(log Y+2+N).

No derivative bound is assumed uniformly across arbitrary test
profiles; this is specific to the fixed dyadic receiver.

## 2. High-order unconditional lattice quadrature

By Poisson summation for the smooth compactly supported h_Y,
and N integrations by parts in the Fourier transform, for N>=2,

    |sum_{m>=1} h_Y(m)-int h_Y|
      <= [2 zeta(N)/(2pi)^N] ||h_Y^(N)||_1
      <= ||h_Y^(N)||_1.

The sum over positive integers equals the full integer sum,
because the support is strictly positive. The coefficient is
at most one since zeta(N)<=2 and 2pi>2. This use of zeta at
real integer arguments requires no zero-location information.

Take l=log Y and N=floor(l/10)>=2. Completing the square,

    5N²-Nl = 5(N-l/10)²-l²/20 <= 5-l²/20.

Also l+2+N<=2(1+l). Thus we obtain the explicit estimate

    |sum h_Y(m)-int h_Y|
       <= 2 exp(5) sqrt(Y) (1+log Y) exp[-(log Y)²/20].

This is an unconditional log-Gaussian quadrature reserve, not
a statement about cancellation of prime or Möbius coefficients.
The underlying integral is exactly

    int h_Y = sqrt(Y)[c log Y+c1],

where c=int exp(s/2) C(s)ds, c1=int s exp(s/2) C(s)ds.
In particular c=U(1/2)^2, the smooth-receiver constant, not
the triangular value of RETURN107-E.

## 3. Near-linear cutoff

Put L=log X>=100 and choose

    D=floor X exp(-4 sqrt(L)).

Then D>=X exp(-4 sqrt(L))/2>=1, and for every d<=D,
Y=X/d satisfies log Y>=4 sqrt(L)>=40. The exact finite
Dirichlet convolution Lambda=mu*log gives

    A_C(X)=sum_d mu(d)/sqrt(d) sum_m h_(X/d)(m).

Define M0(D)=sum_{d<=D} mu(d)/d and
M1(D)=sum_{d<=D} mu(d)log(d)/d, exactly as in RETURN107-E.
The small-divisor contribution is

    sqrt X [(c log X+c1)M0(D)-c M1(D)] + R(X,D).

Use |mu(d)|<=1 and the quadrature estimate termwise. The
square-root prefactor after multiplying by 1/sqrt(d) is
sqrt(X)/d. Since (log Y)^2/20>=4L/5,

    |R(X,D)|
      <= 2 exp(5) sqrt(X) (1+L) X^(-4/5) sum_{d<=D} 1/d
      <= 2 exp(5) X^(-3/10) (1+L)^2.

This is o(1), an improvement in the quadrature component
over the preceding O(log X) reserve. It is unconditional.

The remaining contribution is exactly

    H_C(X,D)=sum_{1<=m<=e² X/D} log(m)/sqrt(m)
        sum_{d>D} mu(d)/sqrt(d) C(log(dm/X)).

Since D>=X exp(-4 sqrt(L))/2, the outer index range is at
most 2e² exp(4 sqrt(L)) = X^o(1). The inner range is still
of macroscopic size for small m. The reduction is NOT a
subpolynomial-time computational algorithm and does not
give any signed estimate on H_C. Indeed trivial absolute
bounds still give O(sqrt X (1+log(e² X/D))²), which here
is O(sqrt X log X), not subpolynomial.

## 4. Coupled signed receiver and its actual remaining obligation

The exact discrepancy is

    A_C(X)-c sqrt X = S_C(X,D)+R(X,D),

    S_C(X,D)=sqrt X [(c log X+c1)M0(D)-c M1(D)-c]+H_C(X,D).

The reserve tends to zero. By RETURN107-B, either eventual
one-sided constant bound on S_C at this cutoff is equivalent
to RH, as is either one-sided O_epsilon(X^epsilon) bound for
every epsilon>0. The implication does not require a uniform
bound across kernels. No estimate of this signed combination
at the required scale is proved here.

The key distinction from RETURN107-E is that there is no
longer a fixed power-length outer sum, and the quadrature
error now vanishes rather than merely being subpolynomial.
Neither fact controls the signed inner Möbius sums or the
truncated Dirichlet moments. Any proposed cancellation method
must still address their combined value. Do not count the
cutoff change as evidence for a global sign theorem.

## 5. Governance and next action

The selected smooth kernel was already zero-detecting and
RH-equivalent, and its explicit derivative growth permits
increasing quadrature order with X. This proves a quantitative
component estimate, not another assertion that smoothing
annihilates off-line modes. The latter would contradict the
nonvanishing multiplier identity of RETURN107-B/C.

Next investigate an unconditional estimate of the coupled
signed receiver using its short outer range and long inner
Möbius windows. In particular do not treat the short outer
range as a small norm: the inner sums carry the unresolved
arithmetic cancellation. If no new sign or cancellation
inequality is identified, report that blocker rather than
repeatedly restating RH-equivalent bounds.

No claim of novelty relative to Poisson quadrature or
classical smooth convolution is made. Internal derivation
only; external mathematical review of the local candidate
remains open. RH remains unresolved.
