# RETURN107-C — Dyadic transport and elementary tent receiver

Status: analytic reduction only; no bound on the actual prime discrepancy,
no new numerical run, and no theorem promotion. Reopens and uses the exact
single-profile construction and one-sided criterion of RETURN107-B.

## 1. Normalize the arithmetic discrepancy before smoothing

On the real log-coordinate line put

    D(ds)=sum_{n>=2} Lambda(n)/sqrt(n) delta_{log n}(ds)
          - exp(s/2) ds.

This is a locally finite signed measure, not a finite total-mass measure.
For any even compactly supported continuous kernel K, set

    E_K(L)=(K*D)(L)
          =sum_n Lambda(n)/sqrt(n) K(log n-L)
             -exp(L/2) integral exp(t/2) K(t)dt.

All sums and integrals in this definition are finite on compact sets
of L. In particular E_{C_b}=A_b-c_b exp(L/2), with C_b the
autocorrelation of the dyadic smooth probability density of half-width b.
The sign is the arithmetic discrepancy, opposite to B in RETURN107-B
up to the explicitly bounded separated Gamma residual.

Write V_b=(1/b) 1_{[-b/2,b/2]} and let u_b denote the dyadic density
of RETURN107-B at scale b. The convolution identity is

    u_b=V_b*u_{b/2}.

The even autocorrelation C_b=u_b*u_b therefore satisfies

    C_b=T_b*C_{b/2},
    T_b=V_b*V_b=(b-|t|)_+/b².

Here T_b is an even nonnegative probability density with support [-b,b].
Its bilateral Laplace transform is

    H_T(z)=[sinh(b z/2)/(b z/2)]².

No scaling identity for D is assumed. Associativity, justified on
compact L sets by local finiteness and compact support, now gives

    E_{C_b}(L)=integral T_b(t) E_{C_{b/2}}(L-t) dt.

The main terms must be included in E before this step: indeed
c_b=H_T(1/2)c_{b/2}, not c_b=c_{b/2}. Iterating,

    C_b=T_b*T_{b/2}*...*T_{b/2^(N-1)}*C_{b/2^N}.

Each transfer preserves eventual upper or lower constant bounds
from the finer receiver to the coarser receiver (shift the starting
point by b). Likewise a one-sided bound M exp(epsilon L), epsilon>0,
transfers with factor H_T(epsilon), not a gain smaller than one.
Since T_b is a nonnegative probability kernel, H_T(epsilon)>=1.
An infinite product of these factors converges, but cannot reduce
a positive exponential growth exponent. The reverse transfer is
not an inequality: differentiating in the sense of distributions,

    E_{C_b}''(L)
      =[E_{C_{b/2}}(L+b)-2E_{C_{b/2}}(L)
           +E_{C_{b/2}}(L-b)]/b².

This is a second difference, not a contraction or a bound on
E_{C_{b/2}}. In particular positivity of the smoothing kernel
does not provide positivity of its inverse.

## 2. Remove the infinite convolution from the arithmetic target

Take b=1. The elementary tent

    T(t)=(1-|t|)_+

alone suffices. Indeed C_1=T*C_{1/2}; if E_T has an eventual
one-sided upper or lower bound, then E_{C_1}=C_{1/2}*E_T
has the same type of bound. C_{1/2} is a probability kernel
supported in [-1,1]. By RETURN107-B either bound implies RH.
The same transfer works for the one-sided exp(epsilon L) allowance,
with finite multiplicative factor H_{C_{1/2}}(epsilon).

For the converse assume RH. Let v=1_{[-1/2,1/2]} (unit integral),
so T is its autocorrelation. It belongs to the logarithmic form
domain on every bounded support: its Fourier transform is the
sinc function, decaying as O(1/|r|). To justify the explicit
formula without assuming its smooth version applies to a jump,
take v_eta=v*phi_eta with a real even nonnegative smooth
probability approximate identity supported in [-eta,eta].
On the real axis |F(v_eta)(r)|<=|F(v)(r)| and pointwise convergence
holds. Under RH the standard zero counting bound N(R)=O(R log R)
implies sum_rho |F(v)(gamma)|²<infty by dyadic summation.
For each fixed L, dominated convergence passes the smooth
spectral cross identity to v, since the translation multiplier
has modulus one. On the physical side v_eta converges in the
logarithmic form norm by dominated convergence on Fourier space,
the prime translations are bounded and finite on the fixed
support of v_eta and its L-translate, and the pole terms converge
by compact-support L2 continuity. Thus

    |W(v,v_L)|<=sum_rho |F(v)(gamma)|²=W(v).

For L>1 the separated Gamma residual is uniformly bounded on
any tail L>=L0>1. The exact cross formula of RETURN106-D extends
to v by the same approximation (or direct separated-support
kernel integration). Hence E_T is bounded under RH. This proves
the equivalence for the elementary tent, with no smooth infinite
convolution required in the tested sum.

Explicitly, put X=exp(L), and c_T=H_T(1/2)=16 sinh(1/4)^2.
Each of the following separately is equivalent to RH:

    sum_{X/e<=n<=eX} Lambda(n)/sqrt(n)
          (1-|log(n/X)|) >= c_T sqrt(X)-O(1);

    sum_{X/e<=n<=eX} Lambda(n)/sqrt(n)
          (1-|log(n/X)|) <= c_T sqrt(X)+O(1).

Endpoint weights are zero. Replace O(1) with the one-sided
O_epsilon(X^epsilon) for every epsilon>0 and the equivalence
still holds. The main term is exactly c_T sqrt(X), not sqrt(X)
or 2 sqrt(X). No sign or numerical bound of this sort has been
proved unconditionally here. The usual RH implication uses the
unconditional zero-counting theorem and the explicit formula;
the sufficient direction uses RETURN107-B's Landau argument.

## 3. Elementary cumulative second-difference receiver

Define the first logarithmic Riesz primitive for real L by

    Q(L)=sum_{log n<=L} Lambda(n)/sqrt(n) (L-log n)
              -4 exp(L/2).

It is continuous and locally integrable; the sum is empty below
log 2. Distributionally Q''=D, since the ramp's second
derivative is a delta and (4 exp(L/2))''=exp(L/2). The exact
finite-difference identity is

    E_{T_b}(L)
       =[Q(L+b)-2Q(L)+Q(L-b)]/b².

This follows directly from the piecewise linear ramp or by
integrating Q'' against T_b. At b=1 the main term on the right is
4 exp(L/2)(exp(1/2)-2+exp(-1/2))
=16 sinh(1/4)^2 exp(L/2), consistent with the literal tent
normalization. No unbounded integration constant is silently
discarded: affine terms of Q cancel exactly in the difference.

Thus the fixed-width arithmetic obligation may also be phrased as
one-sided bounded discrete curvature of Q at step one, or one-sided
subexponential discrete curvature. An estimate of Q itself is
sufficient but need not be necessary for a general function.
This identity does not establish such curvature control.

## 4. Why dyadic smoothing cannot itself suppress an off-line mode

For a hypothetical exponential component exp(zL) of the discrepancy,
smoothing by T_b multiplies it by H_T(z). That multiplier has no
zeros for Re z!=0. This is an exact computation for exponential
functions, not an assumed spectral expansion of D. The infinite
dyadic product has the same nonvanishing property by RETURN107-B.
Consequently the multiscale transport does not annihilate the
obstruction detected by the Laplace pole. It may reduce a coefficient,
but cannot change its exponential growth exponent.

There is also no general inference from local smooth averaging
to a one-sided global bound. For a real oscillatory mode
Re exp((sigma+i tau)L), sigma>0, every one of these averaged
receivers remains an exponentially growing oscillation when its
multiplier is nonzero. This illustrates the limitation of the
operation alone, not a claim that the actual prime measure
contains such a mode. The actual arithmetic mechanism must
exclude the obstruction, not merely smooth it.

## 5. Investigation decision

The dyadic recursion is useful as a transport identity but is not
a closing estimate. It permits the receiver to be compressed from
an infinite convolution to an elementary finite tent weighted
prime sum. The next task is to test a genuine arithmetic inequality
on the cumulative primitive or its finite difference, with all
terms retained at the sqrt(X) scale. Any proposed estimate should
state the exact dependence on the moving center X and whether it
is unconditional. Finite numerical checks cannot settle the
unbounded-X condition. No additional local numerical audit is
indicated by this reduction.

No claim of novelty relative to the classical literature is made.
The fixed-aperture candidate and its completed numerical custody
are unchanged; independent mathematical acceptance remains open.
RH remains unresolved.
