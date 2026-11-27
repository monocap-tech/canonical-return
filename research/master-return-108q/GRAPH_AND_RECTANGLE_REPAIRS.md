# RETURN108-Q — Graph notation and inherited rectangle-budget repair

Status: two failed verification receipts reconciled; two isolated source
repairs submitted for native verification. Certified inventory unchanged
at 14/27. No RH or positivity claim.

## Receipts

Graph run 45924941471, job 140263964507, failed at PacketGraphDensity.lean
lines 19 and 57: the notation for nonnegative reals was interpreted as
an order relation on Type, producing LE Type and OfNat Type 0 errors.
The base GraphDensity module built and standalone audit output appeared,
but the combined strict five-name audit did not execute after the packet
file failed. No promotion from this failed job.

Mellin run 45925084429, job 140264418106, failed during dependency
build, before MellinKernelBank.lean elaboration. The actual blocker
was inherited V10B5c.lean line 210 in angular_rectangle_budget:
ring left a polynomial equality with subtype norm and ambient norm
on opposite sides. Source preflight, Elan, dependency resolution,
mathlib cache and evidence upload passed in both runs. The
Mellin module itself has not yet been tested. Reopened actual logs.

## Isolated repairs

Created branch research/return108q-repair from
4e2c8312ff549b11c589669af3750284418a1bed.
PacketGraphDensity.lean: add open scoped NNReal; commit
ed7d208bc60f754784a4f8e14776f7d36e26119b,
blob 2b5554dd2dc5bd8a5114f286114033bec60f68c7.
V10B5c.lean: replace the ring-only proof at line 210 with explicit
change to the same ambient FullL2 norm in all polynomial terms
followed by ring; commit c889b43c842019238708d78e2a7302ac3fc39774,
blob 2c6c3920b6a1744a156ba425650d467fb4d20f27 (reopened).
This is legitimate because V08I.logCarrier is the kernel submodule
of FullL2 with induced norm. No theorem statement changed.
Neither repair has yet been validated by native compilation.
No original RECON branch modified; this is an isolated patch of an
inherited dependency, not a claim about its current independent head.

## Verification activation

Graph workflow .github/workflows/return108q-graph.yml commit
41c5c68e4cc0c1bbdf14902021772185781819fb;
run 45926487392 in progress at last query.
Mellin workflow .github/workflows/return108q-mellin.yml commit
8e0561ae553a924a9d91c68f38cd992ec607d7c3;
run 45926491992 queued at last query.
Mellin preflight adds the exact repaired V10B5c blob pin and preserves
the certified abstract and carrier source pins and original literal
Mellin/V10B6d pins. Graph preflight pins repaired packet file
and unchanged GraphDensity. Strict named audits 5 and 8 respectively.
The failed receipts and prior certified sources remain immutable.

## Next

Read both actual native results. If inherited RECON dependencies
continue to block the Mellin module, consider splitting out its
finite observation and normalized-pair translation proof from the
sampling-budget imports (the terminal statement takes literal
coordinates and compact signed identity as explicit premises).
No such refactoring performed here. Do not count a source repair
as certification. Full graph-density transport and signed defect
attachment remain unresolved, as does RH.
