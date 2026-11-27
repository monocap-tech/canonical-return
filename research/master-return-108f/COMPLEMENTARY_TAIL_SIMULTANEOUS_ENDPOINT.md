# RETURN108-F — Historical complementary-tail restoration and simultaneous endpoint

Status: historical statement recovery plus internal analytic deduction.
No Lean build or RH claim. Corrects the application boundary of 108-D,
not its fixed-background theorem. The exact native signed-column
identification of 108-E remains conditional.

## 1. Recovered archive and scope correction

Reopened Library "Pasted text.txt" (file search result turn674file0,
read first 1000 lines). It records the Prove2Me mission
"Connes–Weil Green realization: canonical Dirichlet carrier and
actual-zero synthesis", including these named reached milestones:
- ConnesRZ.actual_zero_reciprocal_square_summable: complete
  critical-strip zero subtype, analytic multiplicities;
- ConnesGreen.actual_zero_diagonal_marker_recovery;
- ConnesGreen.actual_zero_half_bound_iff_finite_restoration;
- actual_zero_uniform_marker_recovery.

The archived statement at approximately lines 583–699 fixes P
(all positive), M (selected negative), L=PP*, and B_F (complete
COMPLEMENTARY negative synthesis outside a finite packet F
containing the selected packet). For each epsilon>0 the ordinary
half-bound holds iff for each 0<alpha<1 there is F with
||B_F B_F*||<=alpha epsilon and restored marker >=(1/2-alpha) I.
The next milestone states a single F works across 0<t<=T for a
fixed finite upper aperture T and fixed epsilon,alpha; it is
reflection closed. These are archived theorem statements labeled
PROVED, not newly verified executable proof bodies. The archive
confirms full critical-strip indexing and analytic multiplicities
at the statement level and expressly retains reflection duplication
and pair normalization. It does not expose the imported predicate
bodies or exact signed-column code. The historical Git pin remains
unrecovered. Do not infer literal normalization from prose.

Important correction: B_F decreases as F increases. A persistent
genuine channel is NOT necessarily in B_F, and an exhaustive
sequence F removes every fixed channel eventually. RETURN108-D's
fixed or increasing background obstruction does not apply to
this shrinking-tail operation. Its fixed-B conclusions remain valid.

## 2. Relative marker comparison with no lower bound on P

Let P_t here denote the positive covariance L_t (not the synthesis),
M_t the selected synthesis, and R_F,t=B_F,t B_F,t*.
Use a bounded selfadjoint nonnegative covariance P_t and finite
selected coefficient space. For eta>0 and 0<alpha<1, suppose

    0<=R_F,t<=alpha eta I.

Put Q=P_t+eta I, Q_F=Q-R_F,t. Then

    (1-alpha) Q <= Q_F <= Q,

because alpha eta I<=alpha Q. Both operators are uniformly positive
even if P_t is compact or has kernel. Inversion yields

    Q^-1 <= Q_F^-1 <= (1-alpha)^-1 Q^-1.

Define C_0=(I+M_t* Q^-1 M_t)^-1 and C_F similarly using Q_F.
Compression preserves order. Addition of I and inversion give

    (1-alpha) C_0 <= C_F <= C_0 <= I.

For the lower bound use
I+M*Q_F^-1 M <= (1-alpha)^-1 (I+M*Q^-1 M).
Thus 0<=C_0-C_F<=alpha C_0<=alpha I, so the norm error
is <=alpha. This reproduces the conservative archived recovery
estimate and additionally retains a relative comparison. It
does not give C_0>=I/2. No uniform inverse bound on P_t was used.

## 3. A simultaneous shrinking-tail path is admissible

Take any t_n>0 tending to infinity, eta_n>0 tending to zero,
and alpha_n in (0,1) tending to zero. Apply the archived
uniform recovery statement for upper aperture T=t_n and
parameters eta_n,alpha_n, with a fixed selected packet.
This gives a packet F_n. Alternatively, in the mathematically
reconstructed model the Hilbert-Schmidt column bound of 108-E
ensures ||R_F,t_n|| <= sum outside F of squared negative
column norms ->0 at each t_n. Choose F_n large enough to meet
the relative bound. Enlarging F_n preserves the bound (negative
covariance tail decreases), so one can require nesting, reflection
closure and exhaustion of the countable genuine zero index.
This is a mathematical existence argument, not a numerical
cutoff or a practical complexity bound. Positivity of Q_F
follows with lower bound (1-alpha_n) eta_n I.

Consequently an admissible positive simultaneous path DOES
exist for the complementary tail, unlike for a fixed genuine
background. Along every path satisfying the displayed relative
tail bound, ||C_F-C_0||<=alpha_n ->0. Thus the question is the
ordinary marker on a simultaneous support/regularization path,
not an extra tail sign obstruction.

## 4. Ordinary simultaneous endpoint in the attached model

This section assumes 108-E's native-to-raw signed identification
and ordinary Green energy carrier interpretation. The source
pairings are consistent under nested support windows; native
H0^1 inclusions by zero extension are isometric for the
Dirichlet norm. No uniform native coercivity is needed.

For a fixed finite selected load u define
d_t(u,eta)=<M_t*(P_t+eta I)^-1 M_t u,u>.
For any uniformly positive Q and load b the elementary identity

    <Q^-1 b,b> = sup_h (2 Re <b,h> - <Q h,h>)

follows by completing the square at h=Q^-1 b. Here b=M_t u.
For h in the native carrier the load pairing is the literal
selected sample pairing and the quadratic penalty is
||A_+ h||²+eta ||h||_D². Both are unchanged on embedded smaller
carriers. Therefore d_t(u,eta) is nondecreasing in t for each
fixed eta>0, and nonincreasing in eta for each fixed t.
The ordinary marker C_0(t,eta)=(I+d_t)^-1 is consequently
nonincreasing in t in Loewner order, nondecreasing in eta.
Here d_t denotes also its bounded finite-dimensional
representing operator; no identification across physical
adjoint metrics is presumed without the isometric inclusion.

At a fixed T, the eta-down-to-zero ordinary marker converges
in norm to Gamma_T by the ordinary Picard-to-gap identity
(001F, V05B) and the range-closure equality of 108-E. For every
t_n>=T,

    0<=C_0(t_n,eta_n)<=C_0(T,eta_n).

For a nonempty finite genuine selected negative packet,
108-C gives Gamma_T ->0 in norm as T->infty. Given any
epsilon>0 choose T with ||Gamma_T||<epsilon/2; then choose n
large enough that t_n>=T and
||C_0(T,eta_n)-Gamma_T||<epsilon/2. The order inequality gives

    ||C_0(t_n,eta_n)||<epsilon.

Hence for EVERY simultaneous t_n->infty, eta_n->0 path,

    C_0(t_n,eta_n)->0.

On the admissible restored paths of Section 3,
0<=C_F<=C_0, so C_F->0 as well. The relative error estimate
also shows the two sequences cannot separate in the limit.
This is a squeeze argument, not an unjustified interchange
of support and regularization limits. It is conditional on
existence of the genuine selected packet and the exact native
attachment, and does not assert an actual off-critical zero.
Under RH there is no genuine nonempty negative packet.

## 5. Disposition

Historical complementary-tail restoration is legitimate at each
relative scale and on suitable simultaneous paths. It cannot
produce a different positive global endpoint in the canonical
attached raw/Dirichlet model. A claimed global half-bound on
a genuine packet conflicts with the derived zero endpoint;
a half-bound on a restricted aperture is not contradicted.
The historical endpoint theorem in the archive explicitly
does not assert the required half-bound unconditionally.

The remaining executable-source obligations are unchanged:
exact imported IsTest/IsCriticalZero/zeroMult bodies, signed
duplicate normalization, and source-level identity of the
historical implementation. The recovered mission provides
stronger statement-level provenance, not proof custody.
No workflow activated, no historical code modified, and no
theorem promotion. Next: investigate reflected duplicate
normalization and whether the statement-level source can
be tied to a surviving native declaration or archive.
RH remains unresolved.
