# RETURN107-A — Two-copy transport and the Laplace pole obstruction

Status: analytic reduction, not a new positivity theorem. This note uses
the exact separated-support formula of RETURN106-D. It identifies a
smaller equivalent global target, not an established bound. No numerical
search or assumption of fixed-aperture theorem acceptance is required.

## 1. Two-copy target

Fix a>0 (in particular a=5/4). For real u in C_c^infty((-a,a))
write C(t)=integral u(x)u(x+t) dx, B(L)=W(u,u_L),
u_L(x)=u(x-L). Translation invariance gives

    W(u +/- u_L)=2W(u) +/- 2B(L).

Thus positivity on both signs of two translated copies is equivalent
to |B(L)|<=W(u). The claim to be investigated is only this inequality
for every such u and every L>2a, not positivity on all multi-window
Gram matrices. We show that even the weaker requirement that B(L)
be bounded for each u as L tends to infinity implies RH.

## 2. Transform the exact separated interaction

Let b=2a and choose L0>b. Put

    H(z)=integral_{-b}^b exp(zt) C(t) dt
        =(integral exp(zx)u(x)dx)(integral exp(-zx)u(x)dx).

H is entire; H(1/2)=M_+(u)M_-(u). For Re z>1/2 define
J(z)=integral_{L0}^infty exp(-zL) B(L) dL. The exact formula
of RETURN106-D has pole term exp(L/2) H(1/2),
prime term -sum_n Lambda(n)/sqrt(n) C(log n-L), and residual
Gamma term -R_L(u,u).

For the prime term first integrate over all real L, then remove
L<L0. Absolute convergence of the sum of integrals for Re z>1/2
follows from sum Lambda(n)n^(-Re z-1/2)<infty and compact support
of C. With t=log n-L,

    integral_R exp(-zL) C(log n-L) dL = n^(-z) H(z).

The omitted integrals over L<L0 involve only n<=exp(L0+b);
each has compact L-support. Their sum

    E0(z)=sum_{n<=exp(L0+b)} Lambda(n)/sqrt(n)
                  integral_{-infty}^{L0} exp(-zL) C(log n-L)dL

is entire. Consequently, using the absolutely convergent Euler
product logarithmic derivative in Re(z+1/2)>1,

    J(z)= H(1/2) exp(-(z-1/2)L0)/(z-1/2)
          + H(z) zeta'(z+1/2)/zeta(z+1/2)
          + E0(z) - G0(z),

where G0(z)=integral_{L0}^infty exp(-zL) R_L(u,u)dL.
The explicit Gamma estimate of RETURN106-D makes G0 holomorphic
on Re z> -5/2; local uniform convergence follows on compact
subsets. The bound also avoids any touching-support singularity.

The z=1/2 pole cancels exactly: zeta'/zeta has residue -1 at
s=1, whereas the first term has residue H(1/2).
At z0=rho-1/2 for a zero rho of multiplicity m with Re rho>1/2,
the right side has residue m H(z0). All other terms are
holomorphic there. This is a literal pole obstruction, not a
prime-number-theorem error estimate.

## 3. Consequence of bounded translated interaction

Suppose B(L) is bounded for L>=L0 for every real u supported
in (-a,a). Then its Laplace transform J is holomorphic in
Re z>0 by dominated convergence on compact subsets. The
identity above extends by meromorphic uniqueness from Re z>1/2
to Re z>0. Each alleged zero rho with Re rho>1/2 therefore
requires H(rho-1/2)=0 for every such u.

This is impossible. Given a specific z0, take a real nonnegative
smooth bump of integral one supported in (-epsilon,epsilon),
with epsilon<a. For either sign,

    |integral exp(+/- z0 x) u(x) dx - 1|
       <= sup_{|x|<=epsilon}|exp(+/-z0 x)-1|.

For sufficiently small epsilon the right side is less than 1/2,
so both factors of H(z0) are nonzero. This contradicts
holomorphicity. The standard functional equation symmetry of
nontrivial zeta zeros then excludes Re rho<1/2 as well.

More generally a profile with H(z0) nonzero cannot have
B(L)=O(exp(delta L)) for any 0<=delta<Re z0: the same argument
extends the Laplace transform to Re z>delta. This does not
specify a sign, rate of oscillation or actual off-line zero.

Conversely, under RH the unconditional explicit formula on
smooth compact tests becomes a sum of nonnegative squares on
the diagonal; absolute convergence and Cauchy-Schwarz give
|B(L)|<=W(u). Indeed translation multiplies each real-ordinate
Fourier value by a unit phase. Therefore, for any fixed a>0,
the following statements are equivalent:

    RH;
    every real smooth profile supported in (-a,a) obeys
       W(u +/- u_L)>=0 for both signs and every L>2a;
    for every such profile B(L) is bounded as L tends to infinity.

The second implies the third by the elementary two-copy
identity; the third implies RH by the pole argument.
No assumption of local diagonal positivity is needed for
the third-to-first implication.

## 4. Governance and investigation implications

A universal bounded moving-window discrepancy estimate for
the tested autocorrelations already carries RH-level content.
It must not be slipped in as a harmless technical lemma in
the globalization plan. In particular a summable absolute
cross-block bound is more than sufficient and cannot be
presumed from the fixed-aperture certificates.

The local candidate controls one-copy W(u) on a bounded
aperture; this reduction identifies unbounded translation
as a logically separate obstruction. A single generic bump
might suffice if its bilateral Laplace transform has no zeros
at the relevant zeta zeros, but no such uniform nonvanishing
claim is made here. The argument only needs the whole class
of arbitrarily small bumps within the fixed aperture.

Next global work should identify a genuine arithmetic mechanism
for these tested moving-window discrepancy integrals, or a
different equivalent receiving target. No further local
numerical strengthening can by itself supply that mechanism.

This note uses standard meromorphic continuation, Euler product
and functional equation of zeta, and the explicit formula for
the reverse implication. It is an internal analytic derivation,
not independent review, a Lean proof, or RH acceptance.
