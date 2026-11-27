# RETURN108-W — Combined native build and literal bridge repair

## Reconciled verification
Run 45929439860 job 140278387244 at c741312a475848535ef7ddcf969bfd152388353b. Authenticated closure PASS: 88 modules, 65 IV, 23 RECON/RETURN overlays, 68 distinct external imports, zero dual-lineage required modules. Native targeted dependency build PASS 9010 jobs. Literal bridge audit FAIL at LogCarrierBridge.lean:100:8, a single beta-reduction/rewrite failure in carrierEquiv.left_inv. Seven printed clean axiom sets; carrierEquiv and carrierEquiv_physical inherited sorryAx. Withhold all nine until repaired strict audit. Certified RETURN108 inventory remains 35.

## Independent archive check
Artifact 15168775440 downloaded to local container. Exact archive SHA256 02b62e4f61cbf15928c1aff7e78fe18fd9b2fb99ffc8ec4b2392a5000c6bf66c, matches GitHub digest; ZIP testzip passed. Entries return108v-closure.json, return108v-closure.log, return108v-deps.log, return108v-lean.log. Exact JSON SHA256 ce1f57d6e5d270413167c6aac4e0ed7a9021891c9381b95c69eceeeb960fdb24. Inspected manifest: IV Neutral blob 9c446407e5f9ecfada57d2dee5e1d85b7d8a5753 required_by IV only; original bridge blob 268105f19fd9bb8e1401219268b6bf5d42fef05d; sole cross-edge bridge -> NeutralLogHilbertCarrier. Archive independent integrity confirmed, not independent rebuild. Full fetched job transcript committed here. Exact manifest still in Actions artifact; new workflow prints manifest to log to allow permanent text custody. Do not confuse archive digest with independently reproduced Lean result.

## Repair
After apply iv_physical_injective a, insert:
```lean
    change neutralLogPhysical (toIV a (toRecon a g)).val = neutralLogPhysical g.val
    rw [toIV_physical, toRecon_physical]
```
No theorem statements, hypotheses, supplier files or original lineages modified. Source commit d3067677cad7a0a583ef14b5466018dfad5f0afe, bridge blob 2d45c6dac5353906a0daf2b9d500541be3be3e46. Workflow head 62e893222ed96e3aedae0de5df82b779d5017dda on research/return108s-kernel-core; remote head, source, assembler and workflow reopened. Assembler research/return108w/assemble_closure.py blob 91964b13b7e53e1f22a534b0b7ae46800331c3b1; workflow .github/workflows/return108w-bridge.yml blob e225b6a1c9a504bf16bec8b16a0f982496d23cc0. Pinned repaired RECON source commit and bridge hash; identical IV pin. New run 45934178977 observed. Job state [{"id":140293462771,"status":"in_progress","conclusion":null,"steps":[["Set up job","completed","success"],["Run actions/checkout@v4","completed","success"],["Authenticated minimal closure","completed","success"],["Install Elan","completed","success"],["Resolve pinned dependencies","in_progress",null],["Fetch mathlib cache","pending",null],["Build literal bridge dependencies","pending",null],["Compile and strict named axiom audit","pending",null],["Preserve integration evidence","pending",null],["Post Run actions/checkout@v4","pending",null]]}].

Next: require full nine-target native audit and durable closure receipt, then proceed with squared norm comparison via actual carrier norm identities and continuous equivalence. No norm transport source staged while equivalence is uncertified. RH unresolved.
