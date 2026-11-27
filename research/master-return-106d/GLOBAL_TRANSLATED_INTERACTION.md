# RETURN106-D — Global translated-support interaction: exact receiving object

Status: analytic reduction, not a new positivity certificate. The local
RETURN105-B fixed-aperture candidate remains under literal attachment
review (RETURN106-B2). This note changes the investigation's object
from local diagonal energies to the exact off-diagonal interaction.
It does not launch a runner or assume that the local result is accepted.

## 1. Exact global bilinear form on separated supports

For real compactly supported smooth f,g use the full Guinand–Weil
physical bilinear form, i.e. the form in RETURN098-A with all prime
powers, not the aperture-truncated list of eight. Write
M_+(f)=integral exp(x/2) f(x) dx, M_-(f)=integral exp(-x/2) f(x) dx,
and g_L(y)=g(y-L). Assume f and g are supported in [-a,a],
a=5/4, with L>2a. The gap d=L-2a is strictly positive.
Put C_fg(u)=integral f(x) g(x+u) dx; it vanishes outside [-2a,2a].

The disjoint-support inner product vanishes. Every Gamma increment
therefore contributes only its negative exponential kernel, by
polarization of RETURN098-A. Since y>x for x in supp f and y in
supp g_L, the m=0 term is exactly

    -exp(-L/2) M_+(f) M_-(g).

The two pole terms are

    exp(-L/2) M_+(f) M_-(g) + exp(L/2) M_-(f) M_+(g).

Thus the first pole term cancels the m=0 Gamma term exactly, not
asymptotically. In the prime bilinear sum the shift in the reverse
direction vanishes. The exact result is

    W(f,g_L)
      = exp(L/2) M_-(f) M_+(g)
        - sum_{n>=2} Lambda(n)/sqrt(n) C_fg(log n-L)
        - R_L(f,g),

    R_L(f,g)=integral integral f(x)g(z)
       exp[-(5/2)(L+z-x)]/[1-exp(-2(L+z-x))] dx dz.

For fixed L the prime sum is finite, with exp(L-2a)<=n<=exp(L+2a).
The residual Gamma series is absolutely summable by the positive
gap, justifying the exchanged sums and integrals for signed inputs.
The formula also holds for L2 supported functions with finite
logarithmic form energy, by the same separated-support estimates.
There is no assumption on the zeta zeros.

## 2. Exact prime-discrepancy transport

A change of variables gives

    integral_R exp(s/2) C_fg(s-L) ds
       = exp(L/2) M_-(f) M_+(g).

Let psi(X)=sum_{n<=X} Lambda(n) (Chebyshev function, not digamma),
E(X)=psi(X)-X. With s=log X, the preceding continuous pole term
and the discrete prime term combine to yield

    W(f,g_L)
       = - integral_(0,infty) X^(-1/2)
                         C_fg(log X-L) dE(X) - R_L(f,g).

For smooth compactly supported profiles, integration by parts has
zero boundary terms, and

    W(f,g_L) = integral_R e^(-s/2) E(e^s)
           [C_fg'(s-L) - C_fg(s-L)/2] ds - R_L(f,g).

Define the normalized arithmetic discrepancy
P(s)=e^(-s/2)[psi(e^s)-e^s]. Then the off-diagonal interaction
is a translated compactly supported test against P, minus an
explicit rapidly decaying Gamma kernel. This derivation does not
use a prime number theorem error bound, zero-free region, RH or
any asymptotic assertion about P.

The sign and the factor 1/2 in the differentiated test are fixed
by the derivative of X^(-1/2) C_fg(log X-L), and the upper prime
weight is exactly Lambda(n)/sqrt(n), not 2 Lambda(n)/sqrt(n):
the second prime translation already vanished on disjoint supports.

## 3. Explicit separation-dependent Gamma bound

For d=L-2a>0 the positive scalar kernel
r(u)=exp(-5u/2)/(1-exp(-2u)) is decreasing. Hence

    |R_L(f,g)| <= 2a r(d) ||f||_2 ||g||_2.

This is a rigorous Hilbert-Schmidt upper bound obtained from
the rectangle area (2a)^2. It gives decay in the gap, not a
uniform bound as d approaches zero. Near-touching supports
require a separate bound; do not extrapolate this one to d=0.
The pole and prime terms individually do not satisfy this
decay estimate. Their exact cancellation must be kept before
estimating the interaction.

## 4. Receiving criterion for arbitrary finite separated banks

Let f_i(x)=u_i(x-L_i), with u_i supported in [-a,a] and the
translated intervals pairwise disjoint with positive gaps.
Translation invariance gives W(f_i)=W(u_i).
For any accepted local lower bound c and any symmetric nonnegative
matrix B with zero diagonal satisfying

    |W(f_i,f_j)| <= B_ij ||u_i||_2 ||u_j||_2,

the finite global form obeys

    W(sum_i f_i) >= c sum_i ||u_i||² - v^T B v,

where v_i=||u_i||. A sufficient receiving condition is
||B||_(2->2)<=c, or the stronger row-sum criterion
max_i sum_{j!=i} B_ij<=c. This is a sufficient, not necessary,
criterion and may be prohibitively strong with c=3/10^59.
The exact signed block form can retain cancellations instead.
Even if this separated-bank criterion were proved, a further
partition/density argument would be needed for arbitrary tests,
including touching supports and overlapping windows.

For translates of a single real profile u, the matrix
W(u_{L_i},u_{L_j}) is Toeplitz on a translation lattice
and depends only on L_j-L_i. Full positivity on this span
requires all its finite principal Gram forms to be positive
semidefinite. In particular a necessary two-point condition is

    |W(u,u_L)| <= W(u).

The local aperture certificate gives W(u)>=c||u||², but gives
no upper bound of the required kind for W(u,u_L). Positivity
of each coordinate alone does not yield a positive Toeplitz
operator. A Fourier-symbol argument would require justified
summability or a positive-definite distribution, not merely a
formal infinite series.

## 5. Investigation frontier

The newly specified global arithmetic transport is the family

    integral P(s) [C_fg'(s-L)-C_fg(s-L)/2] ds.

An absolute bound for P is not assumed. The weighted prime
discrepancy on moving windows must be controlled in the actual
test classes, with the right relative energy norm. An arbitrary
fixed finite L census cannot settle the unbounded separation
obligation. The global problem cannot be replaced with another
local Schur optimization; the physical cross object has changed.

Next: finish literal finite-source attachment, then select
a precise signed cross-block or translated-profile receiving
problem. The present note is an analytic identity and
sufficient criterion only. No local theorem acceptance,
arbitrary-support positivity, or RH conclusion.
