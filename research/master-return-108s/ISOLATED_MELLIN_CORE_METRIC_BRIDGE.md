# RETURN108-S — Isolated literal Mellin kernel core and scalar metric bridge

Status: 19 previously certified RETURN108 declarations unchanged.
Mellin bank dependency verification failed again, before compiling its
own source. Isolated smaller core submitted; scalar weight comparison
staged. No RH or positivity claim.

## 1. Actual failed receipt
Run 45926491992 job 140268889502 completed failure.
The corrected V10B5c rectangle budget is no longer the reported
blocking module. Native dependency compilation stopped at V10B5d
line 77, with incompatible complex module instances in
HasDerivAt.const_add/mul_const. The compiler never reached the
eight original MellinKernelBank declarations. Source preflight,
Elan, dependency resolution, cache and evidence preservation passed.
No certification counted from this failed workflow.

## 2. Smaller native receiver
Branch research/return108s-kernel-core from
8e0561ae553a924a9d91c68f38cd992ec607d7c3.
MellinKernelCore.lean commit 3328e22f8fe580a5ebe6447a4aef41589c328a9b,
reopened blob 3364524401994aefe212313bf02ca92cfdb62359.
Imports KernelCompleteness, CarrierTranslation and V10A2
(blob 2dbad9ec81aa2153bacf7039debd51abf69b0102), not V10B6d.
Copies the first three Mellin evaluation theorem bodies from
V10B6b blob 5e8f142d7a86e7b7545606bc9c0e88c49cfe0fff
(carrier_integrand_integrable, carrierMellin_add, carrierMellin_smul).
Corrects the scalar linearity call to the actual argument order.
LiteralCoordinates is the same pointwise proposition as V10B6d,
but a local definition to avoid the upstream imports. The terminal
actual_fullKernel_zero_of_factor receives R,T,K and R.comp T=I+K,
instead of constructing the signed-Gram factor within the module.
Remaining seven candidate statements retain literal sample, carrier
and translation-bank objects. No weakened coordinate premise.
An explicit later signed-Gram adapter is necessary to recover the
original terminal statement; this is a conditional core, not a
proof of existence of bounded sampling or compact defect.
Original eight-theorem candidate preserved. New core includes
eight corresponding core candidates plus three copied helper proofs,
not eleven independent advances over the original eight.

Workflow .github/workflows/return108s-core.yml first committed at
40971796c06ceaf9aad26595bf74cd4dbe75be57. Generated trigger initially
had incorrect branch and file path; corrected in two commits,
final workflow head 4adbbeafa650af98e5a701ac842cb77da9766db5.
Run 45927079259 confirmed in progress. Strict eleven named axiom
audit configured with propext, Classical.choice, Quot.sound.
The original active RECON branch remains untouched.
No claim that core elaborates yet.

## 3. Scalar metric bridge
LogWeightComparison.lean commit
fc55810d684341f388340a40c5a820565731ba90,
reopened blob ded335b61dc283c23ade7357430bc30e275d80ef.
Five new theorem candidates comparing exact expressions
w_I(x)=log(e+|x|), w_R(x)=log(2+|2*pi*x|).
Use explicit additive constants
cRI=log(2+2*pi), cIR=log(e+1/(2*pi)).
For r>=0 and p=2*pi>0,
2+p*r <= (e+r)(2+p), e+r <= (2+p*r)(e+1/p).
Hence w_R<=w_I+cRI, w_I<=w_R+cIR,
and |w_R-w_I|<=max(cRI,cIR), independent of aperture.
Together with w_I>=1 and w_R>=log 2 this also supplies multiplicative
norm comparisons. This scalar variant avoids relying on decimal
bounds for e and pi in the Lean proof. The sharper simple comparison
in 108-R remains valid at paper level.
Five theorem candidates and named print-axioms commands staged;
no workflow activated for this file. These are literal expression
copies, not imports of historical IV091 modules, which are missing
on this branch. The actual cross-lineage same-vector continuous
linear equivalence and compact metric defect remain to be formalized.

## Next
Read the isolated eleven-theorem native audit and repair only
actual reported errors. Then attach the original signed-Gram factor
when its dependency is available. Verify scalar weight comparison
separately and plan the cross-lineage physical and sample integration.
No positivity, RH, or full graph-top declaration promoted.
