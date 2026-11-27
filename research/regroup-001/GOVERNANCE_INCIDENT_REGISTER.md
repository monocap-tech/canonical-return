# Research governance incident register — verification and custody
2026-11-27 Pacific | REGROUP-001 / W02 | Retrospective, append-only companion to R04; not an execution grant, acceptance decision, or amendment of the mathematical program.

## Purpose
Preserve difficulties encountered while conducting reproducible, parallel, modular mathematical research with formal verification. Reproducible mathematics does not by itself imply reproducible research operations. Future investigators may encounter the same classes of failure even with correct proofs and identical sources. Preserve negative operational evidence as well as successful results; avoid treating the observations as inevitable for every implementation.

## Existing evidence and incident taxonomy
R04 already contains a seven-lineage incident synopsis and a proposed governance response:
- [R04 governance protocol, immutable](https://github.com/monocap-tech/canonical-return/blob/505ecf029b849b6d1ca3514252500e3fd43a48c8/research/regroup-001/R04_GOVERNANCE_PROTOCOL.md)
- [W02 instruction precedence, immutable](https://github.com/monocap-tech/canonical-return/blob/e96ebfce615154a5207aafd7b5b6f8e4f99e646e/research/regroup-001/W02_START_HERE.md)

The following classifications consolidate those records and the currently available W02 evidence. They are not an exhaustive forensic audit; report-derived observations are identified as such.

| ID | Lineage / observed problem | Transferable governance lesson | Evidence status |
|---|---|---|---|
| G-01 | RETURN: integration branch triggered ten unrelated historical NF failures; path/ref isolation defects. | Build targets, workflow triggers and mathematical scope must be separately controlled. | R04 report-derived; exact per-run logs not re-audited here. |
| G-02 | IV106: cell-only coercivity argument no longer mathematically needs portions of inherited IV105 response constraints; Lean import/overlay requirements remain. 71 overlay blobs from 35 immutable refs, including 39 generated F08 modules, authenticated. | Logical dependency and compiler dependency are different graphs. A modular theorem need not have a narrow import closure. | [W02 audit preparation](https://github.com/monocap-tech/weil-lab/blob/54b9af8e9b8ea23ec62ce07b3f12f6e9f3b0ea08/verification/iv106/regroup_001/w02_primitive_pole_audit/README.md) read directly. Static preparation only; native audit not executed at that pin. |
| G-03 | NF618: obsolete cumulative ledger digest drove recovery until functional R1 supersession; historical existence unestablished. | Hash receipts are not object custody. Distinguish exact recovery, verified functional replacement and unresolved evidence. | [Frozen custody hold](https://github.com/monocap-tech/weil-lab/blob/afdc8f90715e2da8c703bac1fcf5f9e3aa6e78aa/notes/provisional/NF618CustodyHold.json) read directly; this is a historical state, not a current release decision. |
| G-04 | IV-AM / IV106 / CR-INT: local or expiring archives, hash-only custody, split branches, wrong native commit pointer, binary transfer gaps. | Establish durable transfer capability before computation; authenticate the actual artifact and receiving ref before dependent work. | R04 report-derived. |
| G-05 | CR-ALP: successful workflow with inconclusive terminal test; conditional staging and reporting corrections. | Execution, mathematical validity, custody and adoption are independent status axes. | R04 report-derived. |
| G-06 | RECON / CR-INT: dependent work staged or continued despite failed roots or missing evidence. | Fail closed on affected descendants; independent siblings require documented independence and scope. | R04 report-derived. |
| G-07 | W02: historical instructions, generic continuation permissions, obsolete successor priorities and current grants conflicted. | Explicit instruction precedence and immutable handoffs are part of reproducibility. Repository presence alone is not enforcement. | W02 entrypoint read directly. |
| G-08 | W02: broad historical workflow trigger risk (168 workflows on Weil-Lab main at the historical audit); full isolation not audited. | Preflight actual ref, trigger and path scope before branch or workflow activation. Workflow count alone does not prove that all workflows run. | W02 entrypoint read directly; historical count, not current inventory. |

## Architectural finding: three dependency graphs
1. Mathematical: the declarations, hypotheses and lemmas logically required for the claim.
2. Compiler/build: transitive imports, compiled artifacts, toolchain, cache and actual targets.
3. Provenance/custody: sources, inputs, historical outputs, failure evidence and immutable readback needed to reproduce or audit the work.

Do not equate the three. A cold cache or changed upstream module can amplify compilation over the import closure; a broad workflow trigger can additionally initiate unrelated historical jobs. These are distinct mechanisms, neither of which implies Lean necessarily replays all chronological proof history. G-01 establishes historical cross-workflow amplification by report; G-02 establishes import coupling, not measured rebuild cost.

## Proposed design response (not adopted policy)
- Separate discovery archive (full research history, alternatives and provenance) from certification closure (minimal sufficient theorem dependencies). Preserve both; neither substitutes for the other.
- Define verification scopes V1 local theorem/module, V2 component/interface, V3 cross-component integration, V4 complete clean-room reproduction and named axiom audit. Choose explicitly according to the claim and risk; a scoped build cannot stand in for full certification.
- Use narrow interfaces, isolated pinned source assembly, version-matched compiled dependencies with authenticated origin, and explicit clean-build checkpoints. Treat caches as performance aids, not independent mathematical evidence.
- Measure actual import closure, invalidated/rebuilt modules, workflow fan-out, runner time, cache state and recovery cost. No comparative performance measurements or universal claims are made here.
- Keep a bounded recovery budget and explicit supersession/retirement decisions; do not loop on unchanged evidence. Record first diagnostic, reproduction inputs, attempted remedies and exit conditions.

## Incident entry contract for future investigators
Append a dated record with lineage/operation, expected vs actual behavior, exact source/ref and toolchain, trigger and run/job IDs where available, first diagnostic, affected dependency graph(s), downstream impact, resource cost if measured, source/artifact hashes and readback, attempted remedy, outcome, remaining uncertainty, and transferable prevention rule. Mark observation vs hypothesis vs verified remediation; never invent missing logs or infer mathematical acceptance from a successful CI job.

## Scope and limitations
This initial register is a synthesis of available pinned evidence, not a complete inventory of past or current incidents. In particular, platform access/rate/security errors and suspected lost artifacts reported in conversation require exact diagnostics and receipts before classification as verified incidents. No new mathematics, CI dispatch, workflow change, dependency adoption, or RH/Weil-positivity claim is authorized by this note. No claim of priority or novelty for the general software-engineering concepts.

## Standing owner mandate (2026-11-27)
All present and future investigations and regroups must preserve records and document issues through this register, under the [Standing Incident Reporting Directive](https://github.com/monocap-tech/canonical-return/blob/c8926d30376551cf5f36777828489d7f93b8d355/research/regroup-001/STANDING_INCIDENT_REPORTING_DIRECTIVE.md). The directive establishes reporting, historical backfill, evidence retention, corrections, serialized central indexing, acknowledgement and handoff requirements. Existing source custody and execution restrictions remain intact. Current delivery status: directive committed and remotely read back; separate active thread acknowledgements not yet obtained. This is a governance reporting mandate, not a mathematical progression grant or a technical enforcement claim.
