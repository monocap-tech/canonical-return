# RETURN108-D — Restored-background screening obstruction and regularization boundary

Status: internal analytic deduction from the ordinary raw density theorem
of RETURN108-C; no native Green metric attachment, Lean proof, or RH claim.
Reopened historical monocap-tech/weil main:
- docs/RESTRICTED_CHANNEL_TRANSFER.md blob
  7ee3f7267e5dfc447acf0f9b52e08ceef0487704, WD-B1/B3/B4/B5/B6/B7.
- docs/QUARTET_CHANNEL_RESIDUE_STRUCTURE.md blob
  6a6dad007688e8415f31aa660bb19c5a6644dff5, ZW1-T9/T10 and section 11.
- docs/SUPPORT_FILTRATION_PERSISTENCE.md blob
  513e6fb58b2a892e816e45801484b5df7ea85fb9, terminology of support filtration.
The historical WD-B4 contractive background premise is explicit;
the ZW1 metric separation is explicit. Neither is replaced by
a raw sampling claim. Uses 001G fixed-window positive coercivity
as a paper-level premise for the literal reduced screening map.
The density and negative-witness steps require only RETURN108-C.

## 1. Fixed-window screening norm and exact raw identification

Let X_t be the supported logarithmic carrier with consistent nested
isometric inclusions. Use the actual positive channel A_+, selected
negative packet A_M and remaining negative background A_B, all with
the ordinary multiplicity normalization. B can be finite or infinite,
a fixed set of genuine negative coordinates. The maps are bounded
on every finite X_t. Write S_+=A_+^*, S_M=A_M^*, S_B=A_B^*,
where adjoints use the X_t log metric. Then

    P_t=S_+ S_+^* = A_+^* A_+,
    D_B,t=P_t-S_B S_B^*,
    D_full,t=D_B,t-S_M S_M^*.

For the 001G coercivity premise P_t>=b_t I, b_t>0, P_t is invertible.
The unique reduced solution to S_+ X_B=-S_B is

    X_B,t=-A_+ P_t^(-1) A_B^*.

Indeed S_+ X_B=-P_t P_t^(-1) A_B^*=-S_B,
and its image lies in Ran A_+=(ker S_+)^perp because A_+
has closed range. Any other solution differs by a map into
ker S_+, so this is the reduced solution. Its adjoint on
Ran A_+ sends A_+ f to -A_B f and is zero on the orthogonal
complement. Consequently its norm is precisely

    r_B(t)=||X_B,t||=sup_{f in X_t, f!=0} ||A_B f||/||A_+ f||.

Thus the historical WD-B4 premise in this literal raw realization is
equivalent to r_B(t)<=1, or D_B,t>=0. For a fixed B the ratio
r_B(t) is nondecreasing, since the numerator and denominator
are unchanged on smaller windows. Each fixed-window norm is
finite under 001G; no uniform coercivity constant is assumed.

## 2. Global blow-up for every nonempty genuine background

Choose any genuine negative coordinate j in B. By RETURN108-C,
the full normalized sample vectors of smooth compact tests are
dense in the physical coefficient Hilbert space, with zero
dummy coordinates removed. For each 0<epsilon<1 choose f
with ||A f-e_j^-||<epsilon. It belongs to some finite X_t.
Then

    ||A_+ f||<epsilon,
    ||A_B f||>1-epsilon.

Under fixed-window injectivity f!=0 implies ||A_+ f||>0.
It follows that r_B(t)>(1-epsilon)/epsilon. Given any finite
R>0, choose epsilon<1/(1+R), then choose the finite window
containing this test. Monotonicity gives

    lim_{t->infty} r_B(t)=+infty.

This also holds for a single genuine background coordinate.
It does not rely on a finite-exception assumption, a background
covariance operator on a global physical Hilbert space, or an
exchange of limits. If B is empty or contains only dummy
identically zero coordinates, its analysis is zero and this
claim is inapplicable.

Taking epsilon=1/4 gives an actual conditional existence
statement: there is a compact smooth test f for which

    <D_B,t f,f>=||A_+ f||²-||A_B f||² < -1/2.

No particular off-critical zero or test has been exhibited.
The threshold t is finite but unquantified. By nested consistency,
the same negative form value and log norm persist on every
larger carrier. Therefore once the background fails screening
it cannot later regain screening on a larger window.

WD-B4's residual budget I-X_B X_B^* is nonnegative exactly
when r_B<=1. It therefore necessarily ceases to be nonnegative
on sufficiently large windows conditional on a genuine B.
A positive effective synthesis S_eff=S_+(I-X_B X_B^*)^(1/2)
cannot be defined there by the historical nonnegative square
root. Sequential WD-B6 eliminations cannot fix this: a
legitimate sequence of positive residual covariances would imply
nonnegativity of the aggregate eliminated background form,
which the same argument contradicts if that aggregate is nonempty.

## 3. Regularized restored marker comparison where admissible

For eta>0 and fixed t set

    C_0(t,eta)=(I_M+S_M^*(P_t+eta I)^(-1) S_M)^(-1),
    C_B(t,eta)=(I_M+S_M^*(D_B,t+eta I)^(-1) S_M)^(-1).

These are the ordinary and restored Picard expressions in the
raw metric; no assertion of identification with every native
historical implementation. Require D_B,t+eta I to be uniformly
positive for the second expression and its order arguments.
Since D_B,t+eta I <= P_t+eta I, inversion reverses the order
on uniformly positive operators. Compression by S_M preserves
order; addition of I_M and a second inversion reverse it again:

    0 < C_B(t,eta) <= C_0(t,eta) <= I_M.

For finite selected M, if D_B,t>=0 then as eta decreases to
zero both markers have monotone norm limits (finite-dimensional
positive order boundedness). The inequality passes to the
limits. The ordinary limit is Gamma_M,t by the Picard-to-gap
attachment of 001G/001F. Therefore the legitimate restored
marker is no larger than the ordinary raw marker, not an
improvement over it. This does not identify the effective
analysis range if D_B,t has a kernel or nonclosed range, and
does not assert equality of the two markers.

The sign condition is essential. Mere invertibility of an
indefinite D_B,t+eta I does not justify the displayed positive
order comparison or interpretation as a positive covariance.

## 4. Quantitative obstruction to simultaneous eta->0 and t->infty

For any nonempty genuine B fix the compact smooth witness f
from Section 2. Put

    kappa=(||A_B f||²-||A_+ f||²)/||f||_log² > 0.

For all t containing its support,

    <(D_B,t+eta I) f,f>=(eta-kappa)||f||_log².

Hence a necessary condition for positive restored covariance
is eta>=kappa (strict > kappa for uniform positivity). The
constant is independent of subsequent t, not a uniform
lower bound over all possible B. Thus no path t_n->infty,
eta_n->0 with this fixed genuine B can stay in the positive
restored regime. This rules out an admissible simultaneous
path rather than presuming interchange of two limits.
At any sufficiently large fixed t the positive inner
eta->0 path also fails. There can be indefinite inverses for
some eta below kappa, but they are outside the positivity
premises; a different indefinite operator theory would need
its own exact attachment and theorem.

If B_F is an increasing family of background packets and it
eventually contains the fixed j permanently, use the same
witness from approximating e_j^-. Additional negative
background contributions only decrease D_B,F, so the same
positive lower bound on eta persists for sufficiently large
t and F. This covers monotone inclusion with a persistent
genuine channel, not arbitrary moving backgrounds with no
fixed member.

## 5. Native historical metric remains a separate interface

The historical ZW1-T9 states compact/Hilbert-Schmidt native
Problem-1 Green synthesis and section 11 explicitly refuses
to transfer unweighted sampling coercivity to that metric.
Our blow-up and marker comparison are exact only when the
historical maps are identified with the literal raw adjoints
on the supported logarithmic carrier (or a genuinely
form-preserving isometric/equivalent transport sufficient to
preserve the relevant norm ratios and signs).

If an exact congruence identifies the native background
defect quadratic form with the raw one on a carrier
containing the negative test, the sign obstruction transfers
without needing equality of operator norms. If only a
restricted subspace is attached, one must prove that the
witness is in its image or approximable with form continuity.
Compact synthesis, a bounded finite-head helper estimate,
or an abstract shorted covariance identity does not by itself
supply that attachment. We have not reopened the original
historical executable Green restoration or established a
native form identity; therefore no claim about all restored
historical endpoint geometries follows.

## 6. Investigation disposition

There are two precise ordinary-raw cases:
- Empty genuine background: no restoration; for any genuine
  selected packet the global marker tends to zero by 108-C.
- Nonempty genuine background: its screening norm diverges;
  legitimate positive restoration terminates at a finite
  support scale. In the admissible regime its Picard marker
  is no larger than the ordinary marker, and no positive
  simultaneous regularization path to eta=0,t=infty exists.

No case provides a positive global endpoint independently of
the actual zero-displacement question. The next useful
investigation is the exact native-to-raw form attachment
and its domain/metric permissions, not another unconditional
positive restored endpoint assertion. No historical source
edited, no workflow activated, no independent mathematical
acceptance or novelty claim. RH remains unresolved.
