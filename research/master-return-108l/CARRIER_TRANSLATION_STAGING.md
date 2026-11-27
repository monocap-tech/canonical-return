# RETURN108-L — Isolated kernel verification and literal carrier translation

Status: native verification activated for the five RETURN108-K
abstract candidates; nine additional translation theorems staged, not
compiled. No RH or positivity result. No active RECON lineage modified.

## Custody

Weil-Lab branch research/return108k-abstract-kernel-staging.
Base 486882068fa495d035e9d2450de0e4f215bf6e6c.
KernelCompleteness.lean source commit d83596db4c12ca9826b8b5168692d19a0ddaabc0,
blob 5661627b07544ca200aac74811292e8debb2dff3, reopened.
Workflow .github/workflows/return108l-kernel.yml commit
11dd8cc36828c6fbd9825b975b803441644e1cf3, run 45924177886,
job 140261573094. At last read, checkout, immutable source preflight,
Elan installation succeeded; dependency resolution in progress.
Mathlib cache, native compilation, strict named audit and evidence
upload pending. Allowed axioms propext, Classical.choice, Quot.sound.
The workflow is scoped to the isolated branch and exact kernel file
and workflow paths; subsequent translation staging does not trigger it.

New WeilDefect/MasterReturn/CarrierTranslation.lean source commit
448e47148c4f65580ce20b1a60aff64bf146f096,
blob 19f9587f886e71dc54e0a16d869b8140d0cfbf4a, reopened.
Nine theorem bodies with nine named print-axioms commands; two
noncomputable definitions (shiftPhysical, shiftCarrier).
No new workflow for this file. No local Lean available.

## Literal construction

Inspected V08I, V08J, V08K, V08H, V09B, V07A, V10A, V10B6d
at the pinned branch base; actual logWeight = log(2+|2*pi*xi|).
Inspected Mathlib v4.34.0 Group/Measure, LpSpace/Basic,
LpSeminorm/Basic and Fourier/FourierTransform.
Use translation x->x-s as measure-preserving, Lp.compMeasurePreserving,
and support bound |x|<=|x-s|+|s|. For a+|s|<=b the translated
physical representative is supported in [-b,b] and integrable.
Use the existing Mellin translation identity V07A.mellin_translate
and ae Plancherel/classical compatibility V08H.l2_fourier_integral_ae
with V09B.mellin_imaginary_axis and V08H.angularFourier_scaled.
On the imaginary axis the exponential translation phase has unit
modulus. Thus the literal spectral magnitudes agree ae, the weighted
log energy is integrable, and V08J.liftPhysical gives the actual
translated logCarrier b vector. Exact physical identity, Mellin
covariance and equality of log norms are the intended conclusions.

Staged declarations:
shift_preserving, shiftPhysical_ae, shiftPhysical_support,
shiftPhysical_mellin, shiftPhysical_fourier_norm_ae,
shiftPhysical_log_integrable, shiftCarrier_physical,
shiftCarrier_mellin, shiftCarrier_norm.
They have not been elaborated; API or tactic repairs may be needed.
In particular this is a vector-level shift construction, not yet a
continuous linear map between the two carriers.

## Next

Read actual native kernel audit, repair any failures; separately
verify CarrierTranslation with V09B dependencies. Then attach
linear Mellin observations to Vandermonde receiver and normalized
pair zero transport to the full kernel, using literal V10A fixed
and nonfixed conventions. The full compact-defect receiver remains
conditional on the actual explicit-formula/compact identity supplier.
Green density and full graph transport are not yet native. No
historical ConnesGreen source recovery or RH claim.
