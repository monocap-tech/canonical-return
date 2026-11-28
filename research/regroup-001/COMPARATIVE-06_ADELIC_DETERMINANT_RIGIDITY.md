# COMPARATIVE-06 — Theorem C and adelic determinant rigidity
2026-11-27 Pacific. Owner directed continued Alpöge comparisons only. Source-level analysis, no change to W02 assignments, no external proof adoption, no CI dispatch or native execution. Predecessor canonical-return@0fd79911f55b36dba130a46fcfc411f18470aaee:research/regroup-001/COMPARATIVE-05_SPECTRAL_RANK_SIGNED_DETERMINANT.md.

## Source provenance
Compact rewrite lsernau/grh-lean@72fead102f663d5c7a407964a13dc0cf9d977207, attributed to Alpöge/Claude proof and Sernau compact rewrite per COMPARATIVE-04. Paths below relative to grh-compact:
- src/New/Endgame/LongForms.lean 9eff36cf003ac398fbd24361e80e95a03d0b0e05
- src/New/Heights/Heights.lean a18a99e1f76a6e963113cae8f906ddc99814f518
- src/New/Heights/OneDeterminant.lean 329f004f845b9e2becc896419c8594e82dfef0a6
- src/New/Rows/LiftVolumeHeight.lean 87938141d4103b485fd5a7125721798f65f973f1
- src/New/Rows/HeckeValueHeight.lean f61d2e3c521863d4f814e3338e37c4e2c2fd1f84
- src/New/Selmer/FitTwistsExist.lean da769eba595af2f883eed9b82dbd322a521aea14
- src/New/Rows/RowVolume.lean 925de9c592b40763516e70ddeae48b4ad1a35b88
- src/New/Exchange/CrudeRows.lean e84f4b0e67ed7dbfd727f8ec2dbf05a8d1e5440c (endpoint and relevant excerpt inspected).
These are exact source declarations, not independently accepted transitive proof audits or originality claims.

## 1. Product formula receiver
HeightLe c y is not a naive sum of heights of all y_i. It asserts one common H!=0, archimedean cap alpha and finite-place caps e_l, such that |sigma H|<=exp(alpha) for all complex automorphisms, |iota(H*y_iota)|<=exp(e_l) for every iota:C≅Qbar_l, and alpha+sum_{prime l<M}e_l<=c eventually in M. Thus y_iota can vary with the place, while the common numerator/denominator witness and summed budget bind them. HeightLe.mul and .prod add caps; .pow scales caps. HeightLe.nonneg for constant family 1 derives algebraicity of H from the p-adic uniform embedding bound and uses the number field norm to reduce to the rational logarithmic product formula. It does not assert all family heights are nonnegative. This is an arithmetic, not an Archimedean PSD, constraint.

## 2. Exact Theorem C accounting
Write D=rowDim, H=mainHeight, L=logB*budget=f*a*(log b)^4, Z(v)=normDet(Theta_v)^2. Fit family theta_iota has nonzero Z. Define X_iota=Z(theta_iota)/Z(1). LongForms:
- heightLe_volume_rows gives HeightLe(V+2delta+D*b) X for any independent almost-spanning family with index cap delta and all conjugate untwisted squared covolumes <=exp(V).
- heightLe_volume_ratio_sq gives HeightLe(2D(2L-H)) X^(-2).
- Multiplying square of first with second yields HeightLe(2(V+2delta+D*b)+2D(2L-H)) 1.
- nonneg and division by 2 yield D(H-2L-b)<=V+2delta. Scale inequality b<=L then gives D(H-3L)<=V+2delta.
Both estimates concern reciprocal ratios of the same lift determinant. Their proofs have different inputs but share fit twists and geometric data; not logically independent in an axiomatic sense.

## 3. One determinant: integrality through congruence
Complete n good rows by k=D-n crude rows. crude_rows produces spanning vectors whose twisted lifts are integral after a common nonzero scalar correction gamma(theta)/gamma(1), and whose untwisted conjugate norms are crudely bounded. k<=4D/b^2 bounds cost of completion. one_determinant uses the untwisted Gram determinant as common witness, and local recombination matrices M at j and j∘conj. Integral pairings imply the transformed Gram determinant is integral in the valuation ring. In schematic form det(transformed G)=conj(det M_left)*det G*det M_right, with two local factors separately evaluated using the conjugate embedding. The lower local bound |j(det M)|>=exp(-e_l)|j(gamma_iota)|^k leads to |j(det G_twist)|<=exp(2e_l)|j(gamma_iota)|^(-2k) in compatible notation (actual proof uses conjugate pair of maps). After height for |gamma|^-2 and archimedean cap alpha of completed untwisted Gram, the ratio bound is alpha+2delta+k*g. LongForms charges completion and scalar cost by D*b. This is determinant magnitude/divisibility, not preservation of signed inertia; recombinations can depend on place, and full-rank nonzero determinant is essential.

## 4. Lift volume repayment and scale
RowVolume proves Z(v)=|c(v)|^D product_{j=1..m} lambda_j(v), all nonzero, where lambda_j(v)=normDet(goodLift(v^j,j))^2 and the scalar ratio correction has bounded height under fit twists. This factorization is an actual theorem using common Hecke/Rallis inner-product identities, not a generic rule for determinants of rectangular chains.
HeckeValueHeight.heightLe_inv_ratio gives cap 6B0-dirCost per suitable rank-one pair; heightLe_inv_heckeValueProd_div gives cap 10m*x^999+k(6B0-dirCost) for inverse Hecke-value product ratio. The rank-one proof combines integral test forms, forced divisibility (forcedPart l), local denominator bounds and archimedean size. LiftVolumeHeight uses comparison forms and compressed Gram operator determinants with both archimedean and p-adic bounds to derive per-lift HeightLe(D*k*(logB*B0-dirCost)) for lambda_j(1)/lambda_j(theta_iota). Here B0=budget/(m*k), H=m*k*dirCost. Summing gives D(L-H); doubling and adding scalar cost gives 2D(2L-H) above. Thus repetition scales both main term and budget; saving is not free from repetition. It comes from how the admissible parameter estimates and common height cancellation leave the bound m*remainder <= logB^2*budget (COMPARATIVE-04). The exact forcedPart/localAsk origin remains the next source trace.

## 5. Two distinct twist operations
FitTwistsExist.exists_isFitFamily chooses an m-fold automorphic twist separately for each prime/embedding, controlling first m powers, projection denominators and rank-one Hecke ratio divisibility. It invokes bounded Selmer lengths and budgets summed over finite prime sets. This is not the final quadratic prime twist chosen for a common zero-free right half-strip in OneSidedCriterion/ZeroFreeTwist. The first is a finite-place congruence control operation; the second is an analytic residue noncancellation operation. No equivalence claimed. More detailed Selmer and denominator proofs not audited here.

## 6. Comparison and next question
Weil finite physical Gram/congruence machinery aims at sign/coercivity in a specified real metric; Alpöge uses nonzero determinant magnitudes coupled across all valuations with the unconditional product formula. The generic determinant algebra is shared; arithmetic constraints are not. Theorem C's source reveals a preexisting nonnegativity receiver for a constant family, rather than proving positivity of the target Weil form. The exact arithmetic work is to construct the place-dependent witnesses and force divisibility without excess congruence loss; the next comparative target is forcedPart/localAsk and the conversion of prime data into jet divisibility, including where the Selmer budget enters. No independent proof gate or external proof verification completed. No novelty assessment.
