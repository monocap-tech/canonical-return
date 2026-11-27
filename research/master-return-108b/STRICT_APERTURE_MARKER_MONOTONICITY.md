# RETURN108-B — Strict aperture monotonicity of ordinary selected gap markers

Status: internal analytic extension of MASTER-RECON-001G/H, no native Lean
declaration or workflow, no fixed-aperture acceptance or RH promotion.
The original RECON note already excludes persistence at a neutral
contact. This extension concerns all finite selected loads, not just
half-marker eigenvectors. It is conditional on the literal ordinary raw
model and the fixed-window bounds established at the paper level in 001G.

## 1. Source custody and exact model

Reopened monocap-tech/weil-lab at
206c6934536a1fcfa8ec0c661dcf93eb020e5d04:
research/master-recon/001G-infinite-background.md,
001H-analytic-translation-rigidity.md, 001I-terminal-audit.md;
WeilDefect/MasterRecon/V05B.lean, V06A/B/C/D.lean, V07F.lean.
V05B constructs the ordinary positive Gram compression; V06A/B
provide norm-analytic family and analytic inversion with coercivity
premise; V06C/D propagate vanishing right tails for analytic orbits
with a literal real translation identity; V07F attaches finite
selected forcing. None of these is taken as a proof of global sign.
The staged V10B7 integration at this pin failed (RETURN108-A);
no declaration is promoted based on source inspection.

Let X_t be the supported complex log-form Hilbert space on [-t,t].
A=(A_+,A_-) is the actual reflected, multiplicity-weighted analysis
with unitary sum/difference normalization for nonfixed orbits
and the singleton critical-line convention. Take a nonempty finite
packet Pi consisting ONLY of genuine negative coordinates of
nonfixed orbits delta>0, distinct right-half representatives and
positive analytic multiplicities. Critical-line negative dummy
coordinates must be excluded; they vanish identically.

Set A_Pi=(A_+,A_-,Pi) and K_Pi the corresponding coefficient Hilbert
space. ell:K_Pi->C^Pi extracts the finite selected coordinates;
ell* is isometric insertion. The 001G mathematical hypotheses are
consistent restriction on nested carriers, bounded A on each
finite window and positive channel coercivity

    ||A_+ f||² >= b_t ||f||_log², b_t>0.

The latter is obtained there from full explicit formula, compact
L2 defect and superlinear distinct critical-line supply. We use
these as paper-level premises, not new native certifications.
They imply closed R_t=A_Pi X_t and a unique physical preimage.
Let G_t=R_t^perp and Gamma_t=ell P_Gt ell* on C^Pi.

## 2. Exact projection difference and persistence equivalence

For 0<c<t, R_c subset R_t and G_t subset G_c. Orthogonal projections
onto nested closed subspaces commute and their difference is an
orthogonal projection onto G_c intersect G_t^perp. For any u in C^Pi,

    <(Gamma_c-Gamma_t)u,u>
       = ||(P_Gc-P_Gt) ell* u||² >=0.

Equality holds exactly when v=P_Gc ell*u belongs to G_t.
Let k in X_c be the unique vector with

    A_Pi k=P_Rc ell*u=ell*u-v.

If equality holds, v orthogonal to R_t and hence for every h in X_t,

    <A_Pi k,A_Pi h> = <ell*u,A_Pi h> = <u,ell A_Pi h>.

This is the persistence equation for an arbitrary load u. Unlike
the 001H neutral equation, its left side is the ordinary positive
SELECTED analysis norm, including positive and selected negative
channels. It requires no half-marker condition.

## 3. Analytic translation proves k=0

Choose s=(c+t)/2 and r0=(t-c)/4. For real |a|<r0,
tau_a k lies in X_s and tau_-a h lies in X_t for h in X_s.
The common coefficient Hilbert space K contains both full channels.
Let E be the orthogonal block projection retaining all positive
coordinates and the finite selected negative coordinates.
Let D be the bounded displacement generator, block delta*[0,1;1,0]
with zero on critical singletons, norm <=1/2. Ordinate phases
commute with E and D. For real a the translated positive form is

    p_a(f,h)=<A_Pi tau_-a f,A_Pi tau_-a h>.

It is represented on X_s by

    P(a)=A_s* exp(-aD) E exp(-aD) A_s.

The expression defines a norm-entire operator family for complex
a, using the same two exponentials rather than a complex adjoint.
At zero it is coercive by ||A_Pi f||²>=||A_+ f||²>=b_s||f||²_log.
Therefore the inverse exists and is norm-holomorphic in a disk
of positive radius, by Neumann expansion. Shrink the disk to radius
at most r0.

Use the persistence equation with tau_-a h, and put f_a=tau_a k.
Then p_a(f_a,h)=<u,ell A_Pi tau_-a h>. The right hand side is a
finite combination of raw reflected Mellin samples, with real
translation factors exp((-delta +/- i gamma)a) and their
counterparts. Because the load is finite, these have entire
extensions, with conjugated coefficients extended as power
series rather than conjugating complex a. This gives a norm-entire
forcing vector r(a) in X_s, exactly as V06B/V07F.
Consequently f_a=P(a)^-1 r(a) for real a in the disk.

This produces a genuinely holomorphic X_s-valued orbit whose
physical L2 image agrees with real translates of k.
Apply the 001H or V06D compact-support analytic translation
rigidity: the initially vanishing right tail propagates leftward
via bounded tail restrictions and the identity theorem, so k=0.
This argument needs only local inversion for the positive
selected form, not any invertibility of the signed form.

## 4. A finite selected load cannot annihilate the carrier

From k=0, the persistence equation gives <u,ell A_Pi h>=0
for all h in X_t, in particular all smooth compactly supported
complex h in (-t,t). The negative coordinate at a nonfixed
orbit i is

    sqrt(m_i/2) int h(x) [exp((delta_i+i gamma_i)x)
                          -exp((-delta_i+i gamma_i)x)] dx.

Thus the zero load pairing implies that a finite exponential
polynomial with coefficients proportional to u_i or their
conjugates (depending on inner-product convention) vanishes
as a distribution on (-t,t). It vanishes pointwise there by
continuity and identically by analyticity. The 2|Pi| exponents
are all distinct: distinct right-half representatives have
distinct (delta,gamma), and delta>0 separates right and left
exponents. Evaluate derivatives of orders 0 through 2|Pi|-1 at
zero; the resulting Vandermonde matrix is nonsingular.
Since every m_i>0, all u_i=0.

Together Sections 2-4 prove the strict Loewner inequality

    Gamma_c-Gamma_t > 0   for every 0<c<t.

In finite dimension this also supplies a positive minimum
eigenvalue of the difference for each pair c<t, but no uniform
numerical lower bound on that eigenvalue is derived.

This extends 001H's exclusion of persistent half-marker loads
to all finite selected loads in the ordinary raw model.
If dummy zero coordinates or duplicate linear dependencies were
included, strictness on that redundant space would be false.

## 5. Endpoint and continuity consequences

Fix T>0. By the unique preimage in X_T,
intersection_{0<t<=T} R_t = {0}, since an L2 function supported
in the singleton {0} is zero a.e. The decreasing range
projections tend strongly to zero as t down to zero, by the
nested closed subspace projection theorem. Finite selected
compression converts this to norm convergence:

    Gamma_t -> I as t down to 0.

Similarly the right continuity on every c>0 is exactly 001G:
intersection_{t>c} R_t=R_c. Left continuity also holds by
inward dilation in the log norm. For f in X_c, f_lambda(x)
=f(x/lambda) for 0<lambda<1 lies in X_(lambda c);
its Fourier transform is lambda hat f(lambda xi).
Dilation is strongly continuous in weighted L2 with
w(xi)=log(2+|xi|): for lambda in [1/2,1],
w(eta/lambda)<=w(eta)+log 2<=2 w(eta). The corresponding
dilation operators are uniformly bounded on weighted L2
and converge on compactly supported continuous Fourier-side
functions by dominated convergence, hence on the full space
by density. Thus f_lambda->f in log norm. Bounded A_Pi on
X_c gives R_c=closure(union_{s<c} R_s).
Increasing range projections converge strongly from the left,
and finite compression again gives norm continuity.

The finite ordered eigenvalues of Gamma_t are continuous and
strictly decreasing by the finite-dimensional min-max principle
and the strict Loewner inequality. Their limits at t=0 are 1.
Each eigenvalue can meet 1/2 at most once, and if it meets it
at c it is strictly below 1/2 for every t>c. Multiple eigenvalues
may meet simultaneously. The result supplies no contact
location and no derivative bound.

## 6. Exact remaining global supplier

The monotone norm limit exists as t->infty in finite dimension.
Writing R_infty=closure(union_{t>0} R_t) in the selected
coefficient Hilbert space, the same nested projection theorem gives

    Gamma_infty=ell P_(R_infty^perp) ell*.

Strict monotonicity plus continuity and the initial limit I
show that for this packet the absence of finite half-contacts
is equivalent to Gamma_infty>=I/2. If the limit has eigenvalue
less than 1/2, a corresponding ordered eigenvalue has a finite
unique crossing. Equality at infinity is permitted and cannot
by itself produce a finite crossing.

This isolates the exact global endpoint comparison:
for every genuine finite negative packet Pi,
Gamma_Pi,infty>=I/2. By projection duality (001G), that is
equivalent to nonnegativity of every selected signed form on
every finite window; by convergence of the full negative
channel sum and the explicit formula, universal packet
nonnegativity is equivalent to full Weil positivity on smooth
compact tests. This is RH-level, not a closing premise supplied
by 001G/H. The distinction between ordinary raw and historical
restored markers in 001I remains in force.

Strictness is an activation mechanism, not a global exclusion
mechanism. Neither the existence of the limit nor the
monotonicity controls which side of 1/2 it lies on.
The fixed-aperture numerical candidate does not evaluate
Gamma_infty. A subsequent investigation should examine the
actual orthogonal complement R_infty^perp and whether its
finite compressions admit independent arithmetic or spectral
control. No new runner or theorem promotion. RH unresolved.
