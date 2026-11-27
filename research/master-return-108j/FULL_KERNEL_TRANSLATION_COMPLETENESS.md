# RETURN108-J — Critical-supply bypass for full Green completeness

Status: source dependency audit and internal analytic argument; no fresh
Lean build, no RH proof. Reopened RETURN108-I at
6d059191d0522655112fd5513417208e4e4f5bd7. This is an alternative to
its distinct-critical-zero/Jensen proof, not a proof of critical-line
supply or positive-channel coercivity.

## 1. Exact RECON source and verification audit

Pinned master-recon-v10b20-staging:
486882068fa495d035e9d2450de0e4f215bf6e6c.
- V09C.lean blob 10a8d92d3eb464434166298d127d43f743e99482:
  CriticalSample z means Re z=0, Im z>0, zeta(1/2+z)=0;
  SuperlinearSupply P means for every L there exist T>=1 and
  finite distinct S in the radius-T disk satisfying P with L*T<|S|.
  critical_samples_unique and critical_mellin_injective require
  this premise. Native V09C run 45892279980, job 140157580379,
  completed success including "Compile V09C and inspect axioms".
  This is a conditional receiver, NOT certification of the supply.
- V10B5a.lean bd3f40f61bac56deb59a9f1a828d8b266215f519:
  infinite_log_sampling takes count and weighted budget premises.
- V10B5e.lean c49d79e83f8fdd2d2bf03f46b72673bc55dd31ff:
  bounded_arithmetic_samples attaches both reflected sample banks
  to the literal V09D local count bound.
- V10B7d.lean 8973a57e76f820622bfbaecd0d45482b8a13d07b:
  integer_local_count source proof removes count premise.
- V10B7e.lean 266dca24a8d573b9dbe7ce6dd07da459f7168627:
  arithmetic_analysis, signed_receiver, analysis_unique source proofs.
- V10B6i.lean f47dadb234a1750754d39db1769ebc31574292b1:
  positive_kernel needs SuperlinearSupply CriticalSample;
  positive_coercivity additionally needs signed_defect bound.
- V10B8h.lean: arithmetic_core_diagonal uses LiteralEF;
  arithmetic_operator_identity additionally uses SmoothCore density;
  arithmetic_compact_perturbation produces identity+compact;
  signed_defect gives ||g||² <= ||T_+ g||²+M ||physicalInclusion g||².
- V10B9c.lean 20da13260ec677bbfa95040519cd0f091ba9d66f:
  arithmetic_coercivity explicitly takes OldApproximation,
  VerticalStirling, LiteralEF and SuperlinearSupply.
- V07E.lean: samplePair_translate and sampled_translation are the
  exact real translation law (both channels, not positive alone).
The V10B20 head query returned zero runs; earlier V10B7 staging
integration head had one run, previously recorded failed. The
V10B23 activation ref 920294e467ba6308ae1e21c31ce20b162556caea
has only V10B1/B2/B3a in its tree; the branch suffix is not proof
of containing all later staged source. Do not infer certification
of V10B7d/e from these source inspections. Workflow-list endpoint
was rejected (400); this is not a comprehensive receipt census.

## 2. Full sampling kernel finite dimensional on each fixed window

For aperture b>0 let H_b be the supported log carrier, T_b the
literal bounded full two-channel sampling map. Suppose its
signedGram is I+K_b with K_b compact; equivalently use the
signed_defect estimate below and compact physical inclusion.
This is a premise supplied mathematically by the original full
explicit formula, bounded log-symbol residual, and compact
physical inclusion, not a premise of positivity.

If T_b g=0 then signedGram T_b g=0 because it factors as
T_b* J T_b. Thus g=-K_b g. The unit ball of the closed kernel
ker T_b is relatively compact; it is closed and hence compact.
A normed linear space with compact closed unit ball is finite
dimensional. This proves finite dimensionality without the
positive kernel theorem, density of critical samples, or
nonnegativity of signedGram.

Alternatively, if J_b is compact physical inclusion and
||g||²<=||T_b,+ g||²+M_b||J_b g||² holds on H_b, then on
ker T_b (even on ker T_b,+) the norm is bounded by
sqrt(M_b)||J_b g||. Every sequence in its unit ball has a
subsequence whose physical images are Cauchy; applying the
inequality to differences makes the original subsequence
Cauchy. This also proves finite dimensionality. If M_b=0
the kernel is already trivial.

## 3. Translation excludes the full kernel using one enlargement

Suppose f is a nonzero supported log vector on [-a,a] and
its literal Mellin transform vanishes at all centered nontrivial
zeros. Choose b>a. The continuum of translates tau_s f with
|s|<b-a lies in H_b. Translation invariance of the log norm
and the exact transform law F_s(z)=exp(sz) F(z) show all
these vectors are in ker T_b. This uses both channel samples
being zero, so translation mixing cannot reintroduce a channel.

Any finite family of distinct translates of a nonzero compactly
supported L2 function is linearly independent. Proof: if
sum_j c_j tau_(s_j) f=0 in L2 then Fourier transformation
gives F(iy) sum_j c_j exp(i s_j y)=0 almost everywhere
(up to a common Fourier sign). Since f is L1 by support,
its transform is continuous and not identically zero by
Fourier uniqueness. On a nonempty real interval it is nonzero.
The finite exponential polynomial must vanish there (a.e.,
then everywhere by continuity), hence identically by analyticity.
Taking the first N derivatives at a point gives a nonsingular
Vandermonde system on distinct s_j, so all c_j=0.

Thus ker T_b would contain arbitrarily large independent sets,
contradicting Section 2. Therefore the only supported log
vector whose full actual zero samples vanish is zero.
For T_a injectivity itself, use the normalized pair identities
and positive analytic multiplicities, plus reflection closure, to
recover both actual samples from a zero pair vector. All fixed
points carry their original singleton. No distinct critical-zero
supply is required. One larger aperture suffices; no uniform
compactness or norm estimate as b grows is used.

## 4. Smooth-core variant avoids needing the entire compact identity

If the only available explicit-formula attachment holds on
smooth compact tests, this suffices for the Dirichlet uniqueness
needed in 108-I. Suppose h in H0^1(-a,a) has zero transform
at all actual zeros. Zero extend and convolve with a smooth
compact approximate identity. For sufficiently small scale
choose f=h*phi nonzero by L2 convergence, supported in a slightly
larger interval. Its Mellin transform is the product of the
two transforms by compact-support Fubini, so the zero samples
remain zero. It is smooth compactly supported. Place arbitrarily
many distinct small translates in a still larger fixed aperture b.

On their entire algebraic span E the original smooth explicit
formula gives W(g)=0. The original physical Gårding bound
||g||_log²<= W(g)+M_b||g||_2² follows from the log asymptotic
of the gamma symbol and the finite prime/pole contributions
on b. Hence on E it gives ||g||_log²<=M_b||g||_2².
Compactness of the log-to-L2 inclusion on b and continuity
extend this inequality to closure(E), which must be finite
dimensional by the same Cauchy argument. Translation independence
contradicts this. No smooth-core density in H_b or global signed
operator attachment is needed for this particular variant;
the smooth explicit formula, fixed-window Gårding and compact
physical inclusion suffice.

For completeness the compact inclusion at the paper level is
elementary: truncate Fourier frequency at R. The high-frequency
physical L2 norm is <=||g||_log/sqrt(log(2+R)) under the angular
convention. The low-frequency transform restricted to a bounded
spatial window is Hilbert-Schmidt (square integrable Fourier
kernel on the finite rectangle). Thus physical inclusion
is a norm limit of compact maps. No positivity is used in
the Gårding argument.

This gives an independent classical analytic proof of full
Green completeness without Conrey or a distinct zero count,
using classical literal Weil explicit formula and its Gamma
asymptotic instead. It does not assert those premises have
been newly certified natively. It does not replace the
unconditional unit-height upper count needed by the separate
bounded graph lift.

## 5. Exact graph-top target and native dependency plan

In IV091 at f2e8ec89ae4bf8687a652b299b730d6d1baa43f7 the
target is neutralActualZetaGreenClosedSubmodule a ha = top.
108-I's carrier-density transfer is unchanged: Dirichlet
energy density of actual Green columns follows from all-zero
sampling injectivity via the Green energy pairing, then
continuous dense Dirichlet-to-log inclusion and bounded
full-divisor lift give graph density. The finite packet
source graph lift is the same physical vector by uniqueness.
No new Green synthesis inverse is needed.

Two alternative native dependency chains:
A. V09C.critical_mellin_injective + actual SuperlinearSupply
   + Dirichlet pairing and carrier identification.
B. Full signed compact-defect (or smooth-core Gårding) +
   compact physical inclusion + enlargement/translation transport
   + finite independence of distinct translates.
Both then need the actual supported log carrier identification
between IV091 and RECON, the reflected copy isometry of 108-H,
bounded literal sampling on that carrier, dense energy inclusion,
and graph-top closure transport.
Suggested formal declarations (names only, not compiled):
  fullKernel_finite_of_compactSigned
  finite_translates_linearIndependent
  fullSample_injective_of_enlarged_compact
  dirichletGreen_span_dense_of_fullSample_injective
  actualGreenClosedSubmodule_eq_top_of_density
The first is general operator theory, the second general
compact-support Fourier algebra. Their proof premises can be
explicit theorem parameters, not project axioms. The smooth
variant permits a separate supplier from the literal IV
explicit formula without requiring the entire RECON coercivity
theorem. A native cross-import must be checked against both
pinned dependencies and norm conventions before activation.

Crucial limitation: this only proves FULL analysis injectivity.
If T_+ g=0 but T_- g!=0, translation generically mixes the
negative coordinate into positive; the kernel is not invariant.
The argument cannot replace V09C's critical-supply premise
in V10B6i.positive_coercivity. Neither proves positivity.
No runner activated, no source-level declaration newly
certified, no historical ConnesGreen attachment. Next implement
the abstract finite-kernel and translate-independence receivers
in an isolated native branch or reconcile existing modules
that already supply them. RH remains unresolved.
