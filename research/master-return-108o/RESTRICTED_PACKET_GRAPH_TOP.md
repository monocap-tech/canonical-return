# RETURN108-O — Restricted packet bound gives graph-top and surjectivity

Status: analytic refinement and two new generic Lean theorem candidates.
Five previous abstract kernel declarations certified; carrier and graph
audits running. No native certification for the new two declarations.
No global positivity or RH result.

## 1. Stronger terminal receiver

Let G be the complete source graph, X its physical log carrier,
P:G->X the bounded injective physical projection, D the algebraic
Green packet span, S=closure D. Suppose P(S) is dense in X and
||g||_G <= C||P g||_X for all g in D, with a finite aperture-dependent C.
By continuity the inequality extends to S. The restriction Q=P|S is
anti-Lipschitz: apply the bound to differences in S. Since S is complete
and Q continuous, Q(S) is closed in X. Density implies Q is onto.
For any g in G, choose s in S with P s=P g; injectivity gives s=g.
Therefore S=G and P is surjective simultaneously.
This replaces the separate surjectivity premise in 108-N with a
bound only on the algebraic Green packet span. It does not eliminate
the fixed-window sampling estimate or physical density obligation.

In the literal IV091 graph, reopened at
f2e8ec89ae4bf8687a652b299b730d6d1baa43f7:
ActualZetaSourceGraph.lean blob 5a19c9649315241b9e8cdf0dd3f10821cdd8db27;
ActualZetaGreenGraphPackets.lean blob daf2ad553c5a7ebcaa960b4e321679b8bbe21d5e;
ActualZetaRawSamplesToFullSourceGraph.lean blob
66ac4401448f08af37cde97de012609bd190bec9; use exact graph norm and
injectivity declarations recorded in 108-N. On Green packets, with
the literal normalized orbit sample norm squared bounded by
B_a^2 ||P g||_log^2, the raw graph norm identity gives
||g||_G^2 <= (1+2 B_a^2)||P g||_log^2 (unnormalized source P/N
coordinates have twice the canonical total energy). This estimate
uses 108-H's duplicate and multiplicity normalization and requires
exact carrier/sample attachment. It is not the existing signed
energy identity, which cannot bound the unsigned graph norm alone.
No uniform-in-aperture bound needed. No separate construction of
the full-source lift for every log vector required on this route.
Do not claim graph-top or native cross-lineage identity yet.

## 2. Native candidate and verification custody

Weil-Lab branch research/return108m-kernel-repair:
new WeilDefect/MasterReturn/PacketGraphDensity.lean commit
8678492e6217217e98e5662b8c02ac1cc93ad9b4, reopened blob
40f87862ac3c4fad14f620d8da1eb5dff57f4b51.
Two theorem bodies and named axiom print commands:
graph_top_and_surjective_of_closed_packet_bound;
packet_bound_extends_to_closure.
Uses Mathlib v4.34.0 AntilipschitzWith.of_le_mul_dist and
AntilipschitzWith.isClosed_range; exact APIs inspected.
Not yet elaborated. Potential Lean coercion and rewrite repairs
remain possible.

Independent five-declaration graph workflow
.github/workflows/return108o-graph.yml commit
04479a15e2b55bdd780533d9a5cc77075046d292;
run 45924941471 job 140263964507, in progress.
Preflight authenticates both GraphDensity and PacketGraphDensity
blobs and Lean toolchain; PASS at last read. Dependency resolution
in progress; mathlib cache, build and strict named axiom audit
pending. Allowlist propext, Classical.choice, Quot.sound.
The run verifies only mathlib-level graph receivers, not IV091
or literal arithmetic density.

Carrier run 45924632711 job 140262998127 remains in progress.
Preflight, Elan, dependency resolution, mathlib cache PASS;
Build carrier dependencies running; translation compilation and
nine named axiom audit pending. No carrier error log yet.
Kernel run 45924509279 previously certified 5/5.
Inventory: 27 source theorem candidates total, 5 certified,
22 outstanding. No Mellin bank workflow activated before carrier
verification. Active RECON branch unchanged.

## 3. Next

Read carrier and graph workflow receipts, repair actual errors.
Activate Mellin bank verification after carrier supplier succeeds.
For graph-top, prove dense physical Green packet image and restricted
unsigned sample norm bound with literal IV091/RECON conventions.
This avoids redundant full-graph surjectivity construction but does
not discharge explicit formula, signed compact identity, or the
physical Green completeness supplier at native level. RH unresolved.
