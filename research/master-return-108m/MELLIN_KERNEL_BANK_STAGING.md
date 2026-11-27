# RETURN108-M — Actual Mellin translation banks and kernel audit repair

Status: eight additional theorem candidates staged, not compiled.
The first five-theorem kernel audit failed on two Lean elaboration
errors; corrected source and a dedicated verification workflow committed.
No theorem promotion, no Green graph-top theorem, no RH result.

## 1. Native receipt and repair

Weil-Lab run 45924177886 job 140261573094 completed failure.
Source preflight, Elan, dependency resolution, mathlib cache and
compiler evidence upload passed. Compile/strict named audit failed:
- KernelCompleteness.lean:33: ext already reduced the goal to equality
  of ambient values; redundant apply Subtype.ext failed.
- line 71: Fintype.linearIndependent_iff' in this mathlib pin is the
  linear-combination kernel equality, not the finite sum criterion.
  Use Fintype.linearIndependent_iff instead (literal mathlib v4.34.0
  LinearAlgebra/LinearIndependent/Defs.lean inspected).
No mathematical statements changed. No inference of PASS for dependent
theorems printed after failed elaboration.

Repair branch research/return108m-kernel-repair created from
f17f6fe398b2a102ece5e78a5034b55e2257c22b. Corrected
KernelCompleteness.lean commit 92507c9d5e275afa1f254353bebc6b0487cb12e0,
blob 0fb7619245462b9f778d794596187f4e0d1e1443 (reopened).
New workflow .github/workflows/return108m-kernel.yml commit
38ffa7a96080e8c396eaabdd0283a4673bec444a, pinned to corrected blob.
At first query no workflow run was returned. Separate branch avoids
triggering the old workflow with the wrong source preflight.

## 2. Actual Mellin bank candidate

WeilDefect/MasterReturn/MellinKernelBank.lean source
f17f6fe398b2a102ece5e78a5034b55e2257c22b,
blob 68dca5fdf81597acac1fa9bb85f404fda6ef596a (reopened).
Eight theorem bodies and named print-axioms commands, no static
sorry/admit/axiom tokens. One noncomputable observation definition.
Imports actual V10B6d, CarrierTranslation, KernelCompleteness.
Reused V10B6b.carrierMellin_add/smul, V09B.mellin_entire and
carrier_zero_of_mellin_zero, V10A normalizedPair_fixed/nonfixed,
V10A2.pair_energy_sum, V10B6d.LiteralCoordinates, V10B6c.signedGram.

Declarations:
- observation_eq: literal complex-linear Mellin functional.
- nonzero_mellin_neighborhood: nonzero physical carrier vector gives
  a ball of nonzero evaluations.
- finite_meshes: d=(b-a)/(n+1), e=r/(n+1); all n translations
  fit the same enlarged aperture and all n nodes fit the ball.
  Distinct exp(d(z+j e)) nodes follow from their strictly increasing
  moduli for real d,e>0.
- weighted_shift: literal multiplicity-weighted Mellin covariance.
- normalizedPair_shift_zero: full pair zero implies shifted full
  pair zero, separate critical singleton and nonfixed cases.
  No multiplicity positivity hypothesis needed for weighted zeros.
- shift_mem_fullKernel: literal analysis coordinate custody.
- actual_kernel_banks: construct n independent vectors in the exact
  enlarged sampling kernel, via observation matrix row factors
  F(z+j e) and Vandermonde nodes exp(d(z+j e)).
- actual_fullKernel_zero: for A,T literal analyses at a<b with
  signedGram(T)=I+K, K compact, A g=0 implies g=0.
  Factor using R=(P T)^* P-(N T)^* N, where P,N are original
  coefficient projections. No new positivity or density axiom.

All are source candidates, not native certifications. Mesh proofs,
coercions and other API/tactic details may need repair. CarrierTranslation
remains uncompiled. The new bank module has not been activated; only
the repaired abstract module workflow has been committed. Do not
conflate full kernel with positive-channel kernel or claim unconditional
injectivity without supplying the compact signed identity. The native
Green energy pairing and graph-density transport remain separate.

Next: read the corrected abstract verification, repair as needed;
then activate CarrierTranslation verification and, after it succeeds,
MellinKernelBank verification against its pinned RECON dependency.
No modification to the active RECON lineage or historical Green code.
RH unresolved.
