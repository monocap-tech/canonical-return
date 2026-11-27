# RETURN108-A — Cross-route aperture-uniformity obstruction

Status: internal analytic comparison and exact growth estimate. No new
Lean declaration or runner; no acceptance of the local theorem or RH.
RETURN107-G's divisor-quadrature route has reached an RH-strength
Mertens cancellation requirement. This note audits the corresponding
spectral interface, not a claim that either route solves it.

## 1. Literal RECON sources and verification boundary

Reopened monocap-tech/weil-lab at pinned master-recon-v10b7-staging
commit 206c6934536a1fcfa8ec0c661dcf93eb020e5d04:
- WeilDefect/MasterRecon/V07A.lean: mellin, weightedMellin,
  reflected_mellin_translation, sign tau_a f(x)=f(x-a).
- V07E.lean: samplePair, samplePair_translate, sampled_translation.
- V08I.lean: logCarrier(t) is the closed kernel of supportConstraint,
  with physicalInclusion supported in [-t,t] and logCarrier_energy.
  The norm is the angular log weight log(2+|2pi xi|).
- V10A.lean: Index consists of right-half strip zero representatives,
  displacement delta>=0, ordinate gamma, analytic multiplicity m;
  normalizedPair uses a single weighted coordinate and zero negative
  coordinate at delta=0, two unitary sum/difference coordinates otherwise.
- V10B5a.lean: infinite_log_sampling consumes local counts and finite
  weighted rectangle budgets; the latter are explicit hypotheses.
- V10B6c/d.lean: channel projections, signedGram, LiteralCoordinates,
  signed orbit diagonal and conditional positive coercivity.
- V10B6i.lean: full_strip_receiver, positive_kernel and
  positive_coercivity; coercivity requires critical sample superlinear
  supply and compact-defect inequality.
- V10B7e.lean: staged arithmetic_analysis, signed_receiver and
  analysis_unique for each t>=0, not original physical explicit formula
  or positive coercivity. The integration run 45916353107 is completed
  with conclusion failure at this pin. This note does not infer which
  declarations passed from the overall conclusion or promote V10B7e.

The analytic argument below uses the literal mathematical definitions,
the standard unconditional multiplicity/counting facts and the
classical explicit formula where indicated. It does not rely on
successful compilation of the staged receiver.

## 2. Exact positive sampling growth for a single fixed profile

Use u=u_* of RETURN107-B with b=1, real even smooth nonnegative,
support [-1,1], bilateral transform U(z) entire and nonzero
whenever Re z!=0. Let tau_L u(x)=u(x-L).
For a right-half orbit i with a_i=delta_i+i gamma_i,
delta_i>0, define p_i,q_i by the literal RECON normalizedPair.
Writing w_i=sqrt(m_i) U(a_i), evenness and reality give
sqrt(m_i) U(-delta_i+i gamma_i)=conj(w_i).
By V07A/V07E's mathematical translation law, for any L,

    (p_i(L),q_i(L)) =
      exp(i gamma_i L) [cosh(delta_i L), sinh(delta_i L);
                        sinh(delta_i L), cosh(delta_i L)]
                        (p_i(0),q_i(0)).

Equivalently in the unsymmetrized pair the two samples are
exp(a_i L) w_i and exp((-delta_i+i gamma_i)L) conj(w_i).
The sum/difference normalization is unitary, hence the positive
coefficient norm contribution (sum of squared moduli of both
channels, not the signed form) is exactly

    |p_i(L)|²+|q_i(L)|² = 2 m_i |U(a_i)|² cosh(2 delta_i L).

At delta=0 the normalizedPair convention instead gives exactly
m_i |U(i gamma_i)|², independent of L. No factor two is inserted.

The total positive sample energy

    E(L)=sum_i (|p_i(L)|²+|q_i(L)|²)

is finite for each L. Indeed U has rapid decay on the bounded
real-part strip, by repeated integration by parts on the smooth
compact support; combine with standard multiplicity-weighted
zero counting. For each fixed L the exponential translation
factors are bounded by exp(|L|/2). Set sigma=sup_i delta_i,
0<=sigma<=1/2. The formula gives the unconditional upper bound
E(L)<=exp(2 sigma |L|) E(0).

For every off-critical i, m_i>0 and U(a_i)!=0; consequently
E(L)>=m_i |U(a_i)|² exp(2 delta_i |L|). For any 0<d<sigma
choose i with delta_i>d. These inequalities prove, including
sigma=0,

    lim_{|L|->infty} log(1+E(L))/(2 |L|) = sigma.

This is exact and has no cancellation between different ordinates.
It is conditional in the sense that sigma's value is unknown,
not a detection of actual off-critical zeros.

## 3. Uniformity obstruction on the literal logarithmic carrier

The weighted physical log norm is translation invariant:
multiplication of the Fourier transform by a real-frequency
unit phase does not change its weighted L2 energy. For t>=1
take L=t-1, so tau_L u is supported in [-t,t]. It lifts into
V08I.logCarrier(t): take the Fourier transform multiplied by
sqrt(logWeight), which is square-integrable, and unweight it
in physicalMap. This gives a carrier vector g_t of constant
norm N=||u||_log>0. The support constraint is satisfied because
the recovered physical function is tau_L u. This is a direct
analytic lift, not an asserted compiled Lean lemma.

If T_t is any bounded map with LiteralCoordinates t T_t then

    ||T_t|| >= sqrt(E(t-1))/N.

In particular any asserted support-uniform bound on the
positive full-strip coefficient analysis ||T_t||=O(1) would
imply RH. More generally ||T_t||=O_epsilon(exp(epsilon t))
for every epsilon>0 also implies RH. No physical signed
explicit-formula attachment is required for this direction:
only the literal samples, carrier lift and translation law.

Do not claim a converse for the full family of operator norms
under RH without a separate analysis of growth of sampling
bounds over increasing apertures. In fact critical-line
point evaluation can itself grow with support length.
The single-profile trajectory equivalence is exact:
sup_L E(L)<infty iff RH, and its subexponential growth iff RH.

The aperture parameter t in RECON is a half-width, whereas
RETURN's L is a displacement. The translated one-copy profile
fits in t>=1+|L|; the two-copy sum can be recentered into
t>=1+|L|/2. These are different transports, not identical
aperture parameters.

## 4. Relate to the signed physical obstruction

The signed channel difference for each individual orbit is
invariant under common translation. For a nonfixed orbit this
is 2 Re(w_i conj(v_i)), with v_i=conj(w_i); the factors
exp(delta L) and exp(-delta L) cancel. Thus E(L) may grow
exponentially while the signed diagonal remains constant.
This is exactly why boundedness of the translated signed
diagonal does not bound the positive coefficient norm.

For the cross between u and tau_L u, polarizing the signed
Gram and attaching the original explicit formula on a common
carrier gives the RETURN107-A Laplace obstruction. Under the
real even normalization, its literal full-strip terms are
m_rho U(a_rho)^2 cosh(a_rho L) (or equivalently the sum of
m_rho U(a_rho)^2 exp(a_rho L) by zero reflection).
These are absolutely summable for each fixed L. Do not replace
U(a_rho)^2 by |U(a_rho)|² in this signed receiver. The positive
coefficient norm in Section 2 uses modulus squares instead.

RETURN107-B shows that if sigma>0 the cross is unbounded on
both sides, so some u+tau_L u or u-tau_L u is a genuine
negative physical direction at finite L. By recentering,
it lies in an aperture t>=1+L/2. This is a conditional
first-negative-aperture witness implication. It gives neither
a numerical upper bound on that L nor a proof of a first
crossing exclusion. A finite aperture sign certificate can
coexist with a later crossing unless a global exclusion
mechanism is supplied.

## 5. Cross-route dependency decision

RETURN arithmetic: exact prime-discrepancy identity, centered
Mertens integral and classical unconditional estimates exist;
RH-strength signed cancellation remains unproved.
RECON spectral: fixed-aperture literal sample construction and
signed Gram have staged native interfaces; original-form
attachment, compact-defect and positive supply/coercivity
obligations must be distinguished from global sign exclusion.
The aperture-uniform positive sampling estimate is not a
routine technical extension: the single-profile trajectory
already detects all off-critical displacements without
oscillatory cancellation.
Hybrid physical: RETURN local coercivity is a completed
numerical attachment candidate on [-5/4,5/4], not external
mathematical acceptance. Its mass lower bound is not a
support-uniform log-energy bound. Phase IV/F08 and ALP/INT
fixed-domain transfers cannot be substituted for a global
interaction estimate without matching their domains,
normalizations and quantifiers.

Next: investigate the specific first-crossing exclusion
hypotheses on the signed spectral route. A compact-defect
argument at each t, even if valid, does not alone bound T_t
uniformly or prevent a negative direction at larger t.
Do not count unproved uniformity, an RH-equivalent Mertens
bound, or another receiver identity as a closing supplier.
No new workflow activated and RH remains unresolved.
