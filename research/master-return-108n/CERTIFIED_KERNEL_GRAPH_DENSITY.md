# RETURN108-N — Certified abstract kernel and graph-density transport

Status: five abstract Lean declarations native-certified; nine carrier
translation candidates verifying; eight Mellin-bank and three generic
graph-density declarations staged. RH unresolved.

## 1. Corrected native receipt

Weil-Lab run 45924509279, job 140262610737, completed SUCCESS.
Workflow commit 38ffa7a96080e8c396eaabdd0283a4673bec444a.
Kernel source commit 92507c9d5e275afa1f254353bebc6b0487cb12e0,
blob 0fb7619245462b9f778d794596187f4e0d1e1443.
Reopened job logs, all five exact named PASS lines and
"RETURN108-M strict axiom audit PASS 5/5" present.
Allowlist propext, Classical.choice, Quot.sound.
Artifact 15165574502 return108m-lean, 682 bytes,
reported sha256 96e7d288277e118200bdb6a53ed1a198cc13e76d41ea3e57c4f0e4a10c5602cc.
Not independently downloaded or hashed in this pass.
The previous failed receipt 45924177886 is retained.
New certified inventory 5, staged not yet certified 20.

## 2. Carrier verification

On the same isolated repair branch research/return108m-kernel-repair,
added workflow .github/workflows/return108n-carrier.yml at commit
623b783911ea277453f417dc95b9921e7922e028.
Run 45924632711, job 140262998127. Source preflight passed;
dependency resolution in progress at last read.
Preflight pins CarrierTranslation.lean blob
19f9587f886e71dc54e0a16d869b8140d0cfbf4a and V09B.lean
a3ff3f377235e6e297c2753982a5f5cd951b3fde, Lean v4.34.0.
Build V09B dependencies, compile CarrierTranslation and strict nine
named axiom gates, upload dependency and compiler evidence.
Workflow triggers only its own path; no Mellin-bank runner activated.
No edits to active RECON.

## 3. New generic graph-density receiver

WeilDefect/MasterReturn/GraphDensity.lean source commit
c06c66bd57e1cd8770140fa891c9855224e00df8,
blob 7ac84f8fa2ccb3dc746651d2cb420bb576ac8b1f (reopened).
Three complete theorem candidates:
- dense_image_of_dense_packets: continuous dense-range J:V->X
  takes a dense subset D of V to a dense subset of X.
- closed_graph_eq_univ_of_dense_physical: bounded bijective P:G->X
  between complete normed spaces, closed S in G and dense P(S)
  imply S=G. Proof: Mathlib ContinuousLinearEquiv.ofBijective
  gives continuous inverse, so P(S) is closed as inverse preimage.
- closed_graph_eq_univ_of_packet_transport: combine with inclusion
  J(D) subset P(S). No separate graph convergence estimate needed
  once the physical projection is proved onto.
Inspected Mathlib v4.34.0 Analysis/Normed/Operator/Banach.lean.
No compilation/axiom audit for GraphDensity yet.

## 4. Exact historical graph target and boundaries

Reopened IV091 pin f2e8ec89ae4bf8687a652b299b730d6d1baa43f7:
ActualZetaGreenPacketDensity, ActualZetaGreenClosedSubspace,
ActualZetaGreenGraphClosure, ActualZetaHilbertSourceGraph.
The literal graph projection neutralActualZetaHilbertSourcePhysical
is bounded; injectivity already has a source proof
neutralActualZetaHilbertSourcePhysical_injective. Both graph and
physical log spaces are complete. Surjectivity does NOT follow
from injectivity; requires the bounded actual full-divisor sampling
and same-vector full-source lift on every physical log vector.
The existing neutralActualZetaGreenPacketClosedSpan_eq only identifies
the Green packet closure with the Green closed submodule, not top.
Once actual graph projection surjectivity, Dirichlet packet density,
continuous dense Dirichlet/log inclusion, and literal packet physical
attachment are supplied, the new abstract transport yields graph-top.
The original native mixed and diagonal identities on closure would
then apply to the full graph. No positivity consequence by itself.
No cross-lineage carrier identification has been implemented in this
pass; no source equivalence inferred from similar notation.

Next: read the carrier translation audit and repair its actual
compiler errors if any; then verify the staged Mellin bank against
the now-certified abstract receiver and carrier module. GraphDensity
can be verified separately with mathlib only. The signed compact
defect and exact IV091 graph/RECON attachment remain mathematical
and native supplier obligations.
