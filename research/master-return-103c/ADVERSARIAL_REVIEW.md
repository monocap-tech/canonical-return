# RETURN103-C — Adversarial review of fixed-aperture receiving theorem

The RETURN103-B workflow 45915840440 (job 140234713917) passed all
four literal terminal lineages from 11 authenticated custody archives.
This note reviews the mathematical implication rather than searching
for another favorable matrix. It is an internal review, not an
independent referee report or a native Lean verification.

## 1. No circular RH premise in the spectral identification

Reopened the published article Chirre–Gonçalves, DOI
10.1007/s00209-021-02820-9, Proposition 5, not its RH-dependent
Theorem 1 or Lemma 4. Proposition 5 assumes analyticity in a
closed strip wider than |Im z|<=1/2, decay faster than |z|^-1,
and evenness, and sums over the complex ordinates of all nontrivial
zeros. With real compactly supported smooth f and
F(z)=integral f(x)e^(izx) dx, h(z)=F(z)F(-z) meets these hypotheses
by repeated integration by parts uniformly on bounded strips.
Fourier inversion with the paper's 2pi convention gives
hat h(y)=2pi c(2pi y); hence the -log(pi) term is exactly the
constant in g, the prime factor is -2, and 2h(i/2) is precisely
the pole product. No critical-line assumption is imported.
This is a corollary on smooth tests; the physical form is
defined independently on the larger log form domain.

## 2. Audit of lower comparison and normalization

Reopened RETURN096-A constraint_attachment.py and RETURN096-B
symbol_attachment.py. Both parity sectors have the same retained
frequency set of 81 even positive and five odd positive modes;
the constant or sinh adds one. The complement's grid has exactly
the 16302 finite paid even coordinates, followed by the inclusive
tail beginning 32768. Orthogonality to the five odd modes and
sinh (odd parity) is equivalent to the stated bounded Cauchy
functionals by direct integration. The triangular initial-node
partial fractions and invertible Z preserve the kernel, not its
metric. The finite chart mass uses actual C/T0, not a presumed
orthonormal chart.

Recheck the periodization on [-a,a] with a=5/4: since
|x-y|<=D=5/2 and L=5, the nonzero image terms sum to
2e^(-beta L)cosh(beta(x-y))/(1-e^(-beta L)). Pairing with real
f gives the difference of cosh/sinh squares. For odd q the
sinh( x/2 ) constraint removes both the pole and the beta=1/2
negative image. For m>=1, the elementary bound
||sinh(beta x)||²<=e^(beta D)/(4 beta) gives the image norm
bound in RETURN096-B. This sign cannot be imposed on arbitrary
odd f; it is used only on q in P-perp. The retained H and mixed
t retain the pole and image effects, so no omission is inherited
by their matrix blocks. The prime shifts do not wrap on the
period-five circle since 0<log n<D and L-log n>D.

The two complete interval Fourier bases each have nonconstant
squared norm 5/4. The circle has real-mode norm squared 5/2,
and its coefficients are half the interval coefficients.
Thus the factor 5/8 in the periodic form, split baseline .1
and image payment .001 imply A=.099+w/2 and the signed
coefficients 8(v-d)/pi² inside D. The common physical inverse
factor 4/5 follows exactly. No independent normalization of
the chart or error norm is permissible.

The rational lower symbol verifier is conservative: Machin
alternating remainder, positive atanh log series, integer sqrt
brackets for weights, and monotonic g lower. For its cosine
routine, with input |x|<4 at fixed-point scale S, term truncation
errors in units of 1/S obey e_j <= 1+16 e_(j-1)/
[(2j-1)(2j)], e_0=0. e_1<=1, e_2<=7/3, e_3<3 and the
subsequent multipliers are <=2/7, so every error <3 and their
sum <60. The script's 1000-unit reserve is ample. The omitted
alternating tail after j=20 is bounded by 4^42/42!; cosine
is 1-Lipschitz to cover the directed argument interval. This
justifies the fixed-point phase lower comparison without
assuming a platform cosine is correctly rounded.
The remaining g lower follows from sum/integral variation:
for (t+a)/((t+a)^2+b²), a=1/4, b=r/2, variation <=1/b.
The r lower endpoint can be substituted into
log r-1.85-2/r by monotonicity. The small-frequency positive
digamma sum and integral lower are in the actual verifier.
No RH input in this symbol comparison.

## 3. Domain and source transfer

Real-line Gamma increments and periodized increments are
nonnegative as quadratic forms; the finite partial-sum
difference has uniformly bounded image norm, converging in
operator norm. Monotone convergence gives equivalence of their
finite-energy domains without assuming either is finite at
the start. Bounded prime and pole corrections do not change
that domain. Retained functions' zero extensions have
circle coefficients O(1/|k|) by integration by parts, so
the retained projection maps the log domain into itself.
For each fixed retained p, multiplication of those coefficients
by the log symbol still gives an l2 vector. Thus W(p,q)
extends as a bounded L2 functional of q. Partial Fourier sums
need not themselves satisfy the constraint equations. The
constraint kernel is imposed on the full a, where all rows
are bounded. No pointwise convergence or unjustified operator-
domain approximation is used.

Rechecked RETURN098-C's grouped mixed-row identity. The
differential boundary identity gives the cosine mass 2 sigma
v/(r²+ell), sine mass -2 sigma r b/(r²+ell). The m=0
resonance is handled by the finite M=sinh(a)-a at the retained
sinh mode. The m>=1 denominator never vanishes. Prime overlaps
and pole signs use actual parity; the raw c and d are combined
before chart norms. The raw d is divided by its own pole
denominator, not the trigonometric denominator. This is
important in the moment L(x) and was checked in RETURN102-B.
The k=32768 endpoint is counted in the variable mixed model
and inclusive auxiliary tail, with the same A=87/500.

The compensated odd leading tail is not a null source on the
constraint kernel. Its full rank-one cross terms are preserved;
only the global sinh row is null. RETURN090-B and RETURN091-B
include the signed TT and U crosses, RETURN092-C replaces only
K, RETURN093-B replaces only Hraw, and RETURN103-B now proves
literal preservation. The three error norms use the same
constrained inverse, and the receiving step uses triangle,
not a claim of orthogonality. No model error gets a factor
from a different comparator or a different chart.

## 4. Remaining risk and acceptance distinction

The abstract receiving proof is valid given the certificate
premises. This review found no contradiction in the reviewed
domain, reflection, comparison signs, or 4/5 normalization.
It is not an exhaustive independent verification of the
arithmetic implementation. In particular:
- The original finite mixed-row and feature scalar regeneration
  RETURN097-A/B/C still shares mpmath elementary interval
  operations, whereas RETURN101-A/B/C and RETURN102-B use Arb.
- The large limb products in RETURN087 and RETURN092-B/C rely
  on standard binary64 exact-integer dot products with the
  stated absolute partial-sum bounds. A runtime using reduced
  precision would violate that assumption. A finite platform
  spot test alone is not a mathematical proof of every product.
- RETURN103-B proves literal matrix lineage, but not every
  source enclosure from first principles; its PASS cannot
  replace the analytic attachments.
- Native Lean and independent external mathematical review
  have not been performed.

A targeted next acceptance gate is the Arb regeneration of
the entire remaining finite mixed-row and feature scalar bank
using the already validated physical raw coefficients, with
literal containment in the same archived row/feature enclosures.
This replaces the last identified major shared transcendental
primitive dependency instead of resampling terminal margins.
Separately document the finite dot-product trust base (or
replace it with exact integer products if required by a
stronger acceptance standard). Then produce a single
reproducible proof package with explicit arithmetic trust
assumptions and the integrated analytic derivation. Do not
equate absence of a detected flaw in this internal review
with an independent referee acceptance.

No new runner, unconditional theorem promotion, arbitrary-
support extension or RH claim.
