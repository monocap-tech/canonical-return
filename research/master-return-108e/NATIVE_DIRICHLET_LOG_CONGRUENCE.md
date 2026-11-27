# RETURN108-E — Dirichlet-to-log carrier congruence and native witness transport

Status: internal analytic derivation from the archived canonical Green
definitions, conditional on exact original signed index/column attachment.
No native Lean certification, no independent mathematical acceptance,
no RH promotion, no new runner.

Source custody: reopened Library "Pasted text(4).txt", file
file_00000000cbe082308a226a5bff187123, version 1, 229 lines, the
archived Prove2Me canonical model. It defines SupportedTest, windowL2,
energyVector, energySubspace, Physical, sourceLoad, sourceLift,
sourceEmbed, actualGreenOrdinate/Source, rawColumn and RawSynthesis.
Its imported ConnesRZ.IsCriticalZero and zeroMult bodies are absent.
The archived Git pin f02a526d7c7fb8a43356376a2389d6f8afc08e2a could
not be reopened through the Git tree API (422 invalid commit/tree).
Do not treat the archive's named predicates ColumnRealization,
TestPairing or RawSynthesis as proof certificates; their contents
are specifications. No original executable signed-column or
restoration implementation was recovered in this pass.
Also reopened the project file "Formal congruence transfer.txt":
NF662 describes finite real-entry congruence candidates on a
15-generator block; it is not evidence for a global domain map
or certification of the native signed-column normalization.
The historical public ZW1-T9 metric separation and WD-B4 conditions
remain in force. Uses RETURN108-C/D and 001F/G/H.

## 1. Literal projection pairing

Fix t>0 and use the archived model's window measure on [-t,t].
Write E_t(g)=(g',g/2), K_t=closure(span E_t(SupportedTest_t)).
For this analytic interpretation assume SupportedTest includes
the ordinary complex C_c^infty((-t,t)) tests and is exactly that
smooth compact-support class (or has the same energy closure).
This follows in the conventional mathematical model, but the
imported IsTest body is not present in the excerpt and requires
a literal native check. In particular never apply the totalized
windowL2 to an invalid non-L2 source.

For a valid f in window L2 the literal sourceLift(f) is the
orthogonal projection of (0,2f) onto K_t. Since E_t(g) is in
K_t, in the conjugate-linear-first convention,

    <sourceLift(f),E_t(g)> = <(0,2f),(g',g/2)>
                         = int conjugate(f) g.

Both source and smooth test are L2 on the bounded window.
The native frequency is -i times actualGreenOrdinate, where
actualGreenOrdinate=-i(mirror(rho)-1/2), mirror(rho)=1-conj(rho).
Thus the source frequency equals -(mirror(rho)-1/2)
=conj(rho)-1/2. The literal source is exp((conj(rho)-1/2)x).
For rho in the archived index and m=zeroMult(rho),

    <rawColumn_rho,E_t(g)>
       =sqrt(m) int exp((rho-1/2)x) g(x) dx.

This is a column identity, not a claim about the missing imported
index or multiplicity definitions. It does not require an explicit
Green column ODE identity, a Green inverse, or equality of the
Dirichlet and logarithmic norms.

## 2. Native carrier and dense continuous inclusion

Let V_t=H_0^1((-t,t);C) with norm
||g||_D²=int(|g'|²+|g|²/4). On smooth supported tests E_t is a
linear isometry for this norm. Its extension identifies V_t
isometrically with K_t under the test-class premise above:
the energy norm is equivalent to the usual H1 norm and the
definition of H0^1 is the closure of the smooth compact tests.
The image is complete, hence closed, and contains the dense
generating energy vectors. Zero extension of V_t into H1(R)
is valid by approximation by smooth supported functions;
there are no boundary delta derivatives because of zero trace.

Let X_t be RECON's supported log carrier, using angular Fourier
convention with norm squared

    ||g||_log²=(1/(2pi)) int log(2+|xi|) |hat g(xi)|² dxi.

For r>=0, log(2+r)<=2+r<=3+r² (since 1-r+r²>0).
And 3+r²<=12(r²+1/4). By Plancherel,

    ||g||_log² <=12 ||g||_D².

Therefore zero extension followed by the weighted Fourier
identification of V08I defines a continuous injective map
J_t:K_t->X_t, ||J_t||<=sqrt(12), with the same underlying
physical function. Injectivity also follows directly from the
second energy component g/2. The estimate is uniform in t.

The image is dense in X_t, not asserted closed or onto.
To see this, for any g supported in [-t,t] with finite log norm,
dilate inward g_lambda(x)=g(x/lambda), lambda<1. In Fourier
coordinates its transform is lambda hat g(lambda xi). For
lambda in [1/2,1], log(2+|eta|/lambda)<=2 log(2+|eta|);
the dilation operators are uniformly bounded on weighted L2.
Their strong convergence to the identity follows first on
continuous compactly supported Fourier-side functions by
dominated convergence, then by density. Convolve each inward
dilation with a smooth probability mollifier of sufficiently
small support to stay strictly inside (-t,t). The Fourier
multipliers are bounded by one and tend pointwise to one,
giving log-norm approximation by dominated convergence.
Thus smooth compactly supported functions are dense in X_t.

The two norms are not equivalent. For a nonzero smooth compact
profile phi inside the window, g_n(x)=exp(i n x)phi(x) has
||g_n||_D growing like n||phi||_2 while ||g_n||_log²=O(log(2+n)).
The latter follows from log(2+|eta+n|)<=log(2+n)+log(2+|eta|)
and Plancherel. Hence no bounded inverse of J_t on its image.
A dense inclusion does not contradict native compact synthesis.

## 3. Bounded synthesis exists in the reconstructed literal model

Assume the raw index is the distinct nontrivial zeta zeros,
with positive analytic multiplicities m_rho. This is the
required attachment premise for the archived index, not an
assertion that the missing imported predicate was inspected.
For a=rho-1/2, |Re a|<=1/2. On supported smooth g,
Cauchy-Schwarz gives

    |int exp(ax)g(x)dx| <=sqrt(2t) exp(t/2) ||g||_2.

For |Im a|>=1, integration by parts with no boundary term yields

    |int exp(ax)g(x)dx|
       <= sqrt(2t) exp(t/2) ||g'||_2 / |a|.

Since ||g||_2<=2||g||_D, the norm of the corresponding native
raw column is bounded by C_t sqrt(m_rho)/(1+|Im a|) for some
finite C_t. The projection pairing and density in K_t justify
this bound on the column norm itself. The standard
multiplicity-weighted zero-counting estimate N(R)=O(R log R)
implies sum ||rawColumn_rho||²<infty by dyadic summation.
The column synthesis on l2(Z) is therefore Hilbert-Schmidt:
finite sums converge absolutely in K_t by Cauchy-Schwarz,
and the adjoint has exactly the literal column pairings.
This constructs a compact synthesis without invoking the
archived RawSynthesis predicate's separate equality involving
explicit Green energy; that equality has not been proved here.

On X_t use the bounded full raw analysis of 001G at the
paper level. Its composition with J_t is bounded; equality
with the native raw adjoint holds on the dense smooth energy
vectors by Section 1, hence on K_t. Alternatively the previous
column construction proves bounded native analysis directly.

## 4. Exact signed congruence and range closure

Assume the native signed coordinates are identified by a
coefficient unitary with the canonical multiplicity-weighted
reflected Hadamard pairs, one positive singleton per fixed
point, no redundant negative coordinates (or quotient by
the known null duplicates). Require that this identification
also matches the fixed selected/background coordinate sets.
These are the remaining literal signed-column and index
attachment obligations, not hidden consequences of metric
comparison.

Then for each positive, selected negative and background
channel, the native adjoint analysis N_j equals A_j J_t in
the identified coefficient space. The preceding equalities
on tests extend by boundedness. For the signed covariance
operators (native and log metrics respectively),

    D_native,t = J_t^* D_raw,t J_t.

The same congruence holds for the background defect alone,
and separately for the positive covariance. This is a genuine
congruence, not a similarity or equality of operator norms.
Because J_t has dense image, nonnegativity of either bounded
Hermitian form is equivalent to nonnegativity of the other.
The negative index, defined as the supremum of dimensions
of strictly negative subspaces, is also equal: one direction
follows by injection; the other by approximating a basis of
any finite strictly negative subspace in X_t by vectors in
J_t K_t. Its finite Gram matrix remains strictly negative
under sufficiently small norm perturbations by boundedness
of D_raw,t; strict negativity also guarantees independence.

Similarly closure(Ran N_Pi)=closure(Ran A_Pi), since J_t
has dense image and A_Pi is bounded. Therefore the ordinary
gap projection and finite marker coincide under the stated
coefficient identification at each fixed t. No assertion of
closed native analysis range follows. Under 001G closed raw
range, its closure is R_t; a compact native operator may
have a dense proper image in that closed raw range.
The ordinary native Picard endpoint from 001F is the
compression of ker S, which is (Ran S*)^perp; it thus equals
the raw ordinary marker under this attachment, even without
native coercivity or a Green inverse.

## 5. Transport the conditional negative witness

RETURN108-C produces a smooth compactly supported complex
test g with ||Ag-e_j^-||<1/4 conditional on the existence of
a genuine negative coordinate j. It is a finite combination
of the two dyadic translated profiles. Its support is inside
some open bounded interval, so E_t(g) is an actual native
energy vector for all sufficiently large t under the
test-class premise. It is not merely a limit vector in a
larger log carrier. For any background B containing j,

    <D_native,B,t E_t(g),E_t(g)> < -1/2.

This transfers the conditional sign obstruction despite the
stronger Dirichlet norm. The norm of E_t(g) is fixed once the
whole support is inside the window. Hence

    kappa_D = (||A_B g||²-||A_+ g||²)/||g||_D² > 0

is independent of larger t, and a positive native restored
covariance D_native,B,t+eta I requires eta>=kappa_D.
The conclusion does not need a uniformly positive native
positive covariance P; it is a quadratic witness argument.
No simultaneous path t->infty,eta->0 can satisfy positivity
for a fixed genuine B in this exactly attached model. Mere
indefinite invertibility remains outside this claim.

## 6. Remaining literal custody and decision

The carrier metric is not an escape in the canonical
projection realization: the actual smooth witnesses are
already native test vectors, and the stronger norm does
not change the sign. What is left is NOT a global metric
equivalence theorem. It is direct source checking of:
1. imported IsTest and IsCriticalZero/zeroMult definitions;
2. actual native signed-column normalization, duplicate
   quotient and selected/background packet identification;
3. actual executable restoration path, including whether
   it uses this covariance, sign premise and regularizer.
The old source Git pin returned 422; the recovered archive
is a specification, not full executable proof custody.
The NF662 finite congruence note is a separate fixed-block
investigation, not a replacement for these missing sources.
No new runner or historical edit. No theorem promotion or
claim of actual off-critical zeros. RH remains unresolved.
