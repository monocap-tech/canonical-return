# RETURN108-T — Independent scalar audit and literal carrier bridge staging

Status: 19 certified RETURN108 declarations unchanged. The isolated
Mellin core and scalar comparison are both verifying. A cross-lineage
same-physical-vector linear equivalence is staged, but not compiled.
No global positivity or RH result.

## 1. Verification custody

Mellin core run 45927079259 job 140270815878 remains in progress:
immutable source preflight, Elan, dependency resolution, mathlib cache
PASS; build of certified suppliers and V10A2 dependencies running.
No compiler outcome on the eleven candidate declarations yet.

New workflow .github/workflows/return108t-weight.yml at
5660a0bd091579f4ca3eeff959bbadbf4135325f, branch
research/return108s-kernel-core, run 45927305613 job 140271541361.
Pins LogWeightComparison.lean blob ded335b61dc283c23ade7357430bc30e275d80ef
and Lean v4.34.0, builds standalone scalar source against mathlib and
strictly audits all five named theorems, allowlist propext,
Classical.choice, Quot.sound. Source preflight and Elan PASS;
dependency resolution running at last read. No certification yet.

## 2. Exact source-tree integration census

Compared untruncated recursive trees from IV091 pin
f2e8ec89ae4bf8687a652b299b730d6d1baa43f7 (8851 entries) and RECON
return108t workflow head 5660a0bd091579f4ca3eeff959bbadbf4135325f
(827 entries). Counts of .lean blobs:
IV091 620; RECON 195; shared paths 44; identical shared blobs 42;
different shared blobs 2; IV091-only 576; RECON-only 151.
Both lean-toolchain v4.34.0 and identical lakefile.toml blob
e567c529fc39eff8a175b3361a2831763e5dcc2c, mathlib v4.34.0.
Different shared paths:
- WeilDefect.lean: IV 3941843128d7f51a0543366d65b6f20dd5f7fbc6,
  RECON d769b418e3d8215a172398efb2831147c0bfd893 (umbrella imports).
- WeilDefect/Morphology/Neutral.lean:
  IV 9c446407e5f9ecfada57d2dee5e1d85b7d8a5753,
  RECON 5f20024f9d8d6567584f838478dca071848a7618.
  This second conflict is substantive, not just an umbrella import:
  IV091 NeutralArithmeticMorphology has retained densityIntegrable
  and logEnergyIntegrable witnesses absent from the RECON version.
  Reopened both exact files. Do not assume compatible structures.
An integration test should start from IV091 (preserve its actual graph
and source definitions), selectively import RECON and RETURN modules,
authenticate the transitive source closure, reject conflicting
dependency paths unless explicitly reconciled, and run targeted
lake build, not the aggregate root. This is a proposed integration
policy; no overlay or combined build executed. Counts refer to the
specified immutable heads, not the current total project inventory.

## 3. Literal bridge candidate

New WeilDefect/MasterReturn/LogCarrierBridge.lean commit
8866d05c4d795bf153027b0b076a323262b0c71a,
reopened blob 268105f19fd9bb8e1401219268b6bf5d42fef05d.
Eight theorem bodies and three definitions (two lifts and carrierEquiv);
nine named print-axioms commands including the equivalence definition.
The source imports both literal V08J and NeutralLogHilbertCarrier,
the latter absent from this branch. Integration-only; no runner.
Reopened IV091 NeutralCanonicalFormDomain, NeutralLogFormGraph,
NeutralLogHilbertCarrier and RECON V08H, V08J. Use the same Mathlib
physical L2 Fourier representative, not identification of weighted
coordinates. From the additive bounds of 108-S and positive lower
weights derive multiplicative bounds using
C_R=1+|cRI|, C_I=1+|cIR|/log 2. Domination gives energy integrability
in both directions for the identical L2 vector.
toRecon uses V08J.liftPhysical with actual IV091 support and
energy; toIV uses neutralCanonicalToLogHilbert with actual RECON
support and converted energy. Prove both exact physical identities,
then mutual inverse and complex linearity from injectivity.
carrierEquiv is a candidate LinearEquiv, NOT a certified continuous
linear equivalence; norm bounds/continuity and exact sample
attachment remain to be compiled. The mathematical transport is
the one derived in 108-R. No isometry, graph-top, compact defect
or positivity claimed.

## 4. Next
Read both workflow results. Repair actual native errors. Prepare an
isolated combined source closure based on IV091 with explicit conflict
custody; only then compile the literal bridge, its norm bounds and the
actual full-divisor sampling normalization attachment. No original
active branch edited. RH unresolved.
