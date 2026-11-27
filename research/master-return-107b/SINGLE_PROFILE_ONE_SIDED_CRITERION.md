# RETURN107-B — Single-profile and one-sided prime discrepancy criterion

Status: internal analytic reduction; not an arithmetic estimate, independent
referee acceptance, native Lean proof, or RH resolution. Uses the exact
translated interaction of RETURN106-D and Laplace identity of RETURN107-A.
No assumption of accepted fixed-aperture positivity and no new runner.

## 1. Fix an explicit smooth zero-detecting profile

Choose b=1 (or any 0<b<5/4), a_j=b/2^j, j>=1. Let v_j be the
uniform probability density on [-a_j,a_j]. On a product probability
space take independent X_j with these laws. Since sum a_j=b,
sum X_j converges absolutely for every realization. Its law mu is
supported on [-b,b], is even, and has total mass one. The bilateral
Laplace transform is

    U(z)=int exp(zx) dmu(x)
        =prod_{j>=1} sinh(a_j z)/(a_j z),

with the factor at zero defined to be one. Equality follows by
dominated convergence for the finite partial sums, uniformly on
compact z-sets. Alternatively the product converges locally uniformly
since each factor differs from one by O_K(a_j²) on every compact K.

Each factor is nonzero off the imaginary axis. For any fixed such z,
the tail deviations from one have summable absolute values; the
convergent infinite tail product is nonzero (use the power series
log(1+w) for |w|<=1/2). Finite initial factors are nonzero too. Hence
U(z) is nonzero for Re z!=0. This requires no assertion about zeta zeros.

To see that mu has a C-infinity density, for real r each characteristic
factor |sin(a_j r)/(a_j r)|<=min(1,1/(a_j|r|)). For every N>=1,

    |U(ir)| <= (prod_{j=1}^N 1/a_j) |r|^-N,  r!=0.

Thus U(ir) and all r^k U(ir) are integrable (use boundedness near 0).
Fourier inversion for finite measures gives the continuous density
u_* and differentiation under the integral gives all derivatives.
Since the measure is supported on [-b,b], this smooth density
vanishes outside [-b,b]; it is nonnegative and even. For b=1 it
lies in C_c^infty((-5/4,5/4)). Define

    C_*(t)=int u_*(x) u_*(x+t) dx,
    H_*(z)=U(z)U(-z)=U(z)^2.

Then C_* is nonnegative, smooth, supported in [-2b,2b], has integral
one, and H_*(z) is nonzero for Re z!=0.

## 2. One profile suffices

Set B(L)=W(u_*,(u_*)_L), L0>2b. RETURN107-A gives, for Re z>1/2,

    J(z)=int_{L0}^infty exp(-zL) B(L)dL
        = H_*(1/2) exp(-(z-1/2)L0)/(z-1/2)
          + H_*(z) zeta'(z+1/2)/zeta(z+1/2)
          + E0(z)-G0(z).

The first and zeta-pole terms cancel at z=1/2; E0 is entire,
G0 holomorphic for Re z>-5/2. At each nontrivial zero rho with
Re rho>1/2 the residue at z=rho-1/2 is m_rho H_*(rho-1/2),
which is nonzero for this single fixed test. Thus bounded B
would exclude every such zero by Laplace holomorphicity in
Re z>0. Under RH the explicit formula gives |B(L)|<=W(u_*)
as in RETURN107-A. This proves the single-profile boundedness
criterion. No local certificate is used in this equivalence.

## 3. Landau positivity obstruction: one-sided boundedness suffices

Assume either B(L)<=M or B(L)>=-M eventually, taking M>=0.
Increase L0 so that the chosen inequality holds for all L>=L0.
Then F(L)=M-B(L) or M+B(L), respectively, is nonnegative there.
Its Laplace transform is

    K(z)=M exp(-z L0)/z - J(z), or
    K(z)=M exp(-z L0)/z + J(z).

It converges absolutely for real z>1/2 by the same majorants as
RETURN107-A. Let alpha be the infimum of real convergence points,
possibly -infinity. For every Re z>alpha the integral defines a
holomorphic function by local uniform domination (polynomial
factors arising from derivatives are absorbed by any smaller
real convergence point).

If there exists a zero rho with Re rho>1/2, the nonzero pole in
the continued K requires alpha>=Re rho-1/2>0. To justify this,
if alpha<Re rho-1/2, analytic uniqueness on Re z>alpha would
make the pole removable, a contradiction. In particular alpha
is then a finite positive real number, <=1/2.

But the meromorphic expression is holomorphic in a neighborhood
of the real point alpha: zeta has no real zeros for s>1/2.
For 1/2<s<1 this follows from the alternating eta series,
eta(s)>0 and zeta(s)=eta(s)/(1-2^(1-s)); for s>1 use its
Euler product. At s=1 the previously exhibited cancellation
removes the pole of J. Also alpha>0 so the constant transform
has no pole there. Entire E0 and the decaying Gamma transform
are holomorphic there.

Here is the needed Landau lemma, with proof rather than an
assumed growth theorem. A Laplace transform of a nonnegative
locally integrable F on [L0,infty), L0>0, with finite real
abscissa of convergence alpha, cannot admit a holomorphic
continuation to a disk about alpha. Suppose the disk radius is
r>0. Choose sigma=alpha+r/4. The Taylor expansion of K about
sigma converges at sigma-h, with h=r/2, because the disk
centered at sigma of radius 3r/4 is inside the original disk.
For every n,

    (-1)^n K^(n)(sigma)=int L^n exp(-sigma L) F(L)dL>=0.

All derivatives are justified by convergence at a point between
alpha and sigma. The series at sigma-h consists of nonnegative
terms. By Tonelli,

    sum_{n>=0} h^n (-1)^n K^(n)(sigma)/n!
       = int exp(-(sigma-h)L) F(L)dL.

The left side is finite by holomorphicity; the right side would
therefore converge at sigma-h<alpha, contradicting the definition
of alpha. This proves the lemma and the contradiction.

It follows that either eventual one-sided constant bound excludes
all off-critical right-half zeros; symmetry then implies RH.
Conversely RH gives both bounds. Equivalently, if RH is false,
this explicit single B(L) is unbounded above and below on every
tail. This is conditional reasoning, not an assertion of actual
off-line zeros or a computed sign failure.

The same argument permits a weaker growth requirement. If for
every epsilon>0 there is a constant M_epsilon such that eventually
B(L)<=M_epsilon exp(epsilon L) (or the reversed bound), then use
F=M_epsilon exp(epsilon L) - B(L) (or plus B) and the constant
term transform with pole at epsilon. Any alleged rho with
Re rho-1/2>epsilon forces alpha>epsilon, and the argument above
is unchanged. Thus one-sided subexponential growth also suffices.

## 4. Exact arithmetic inequality on one moving window

For b=1 define c=H_*(1/2)=U(1/2)^2>0 and

    A(L)=sum_{exp(L-2)<=n<=exp(L+2)} Lambda(n)/sqrt(n)
                                     C_*(log n-L).

Every summand is nonnegative. The exact separated interaction is

    B(L)=c exp(L/2)-A(L)-R_L,

with 0<=R_L<=2b exp[-5(L-2b)/2]/(1-exp[-2(L-2b)])
                    * ||u_*||_2² for L>2b.
The inequality follows from positivity of u_* and the positive
Gamma kernel. For all L>=L0>2b this remainder is uniformly bounded.
Accordingly each of the following is individually equivalent to RH:

    A(L)>=c exp(L/2)-O(1)  as L->infty;
    A(L)<=c exp(L/2)+O(1)  as L->infty.

One-sided O(exp(epsilon L)) errors for every epsilon>0 are also
individually equivalent. The constants may depend on this fixed
profile and epsilon. In the original X=e^L coordinate the error
is O(X^epsilon) in the weighted sum, and the main term is c sqrt(X).
No uniform in-profile control, signed block estimate, or multi-window
matrix positivity is required by this reduction.

This does NOT supply either arithmetic inequality. The positivity of
the individual summands alone controls neither side to the stated
accuracy. An ordinary relative o(1) prime number theorem bound
on A(L)/(c exp(L/2)) does not suffice: its absolute error can grow
like exp(L/2)o(1). A power-saving estimate with arbitrary epsilon
is exactly the new global obligation, not a free assumption.

## 5. Next investigation

The dyadic convolution gives explicit approximating compact spline
kernels and a functional equation U(z)=sinh(bz/2)/(bz/2)*U(z/2).
It provides a fixed arithmetic target and a possible computational
test family; it does not imply a recurrence for the prime discrepancy
without further arithmetic input. A finite separation census cannot
establish the asymptotic inequality. Seek a genuinely unconditional
arithmetic mechanism for one side of this weighted moving-window
estimate, or a logically weaker but still sufficient condition.
No new workflow was activated and no local or global theorem
promotion is made.
