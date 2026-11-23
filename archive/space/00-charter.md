# Charter — Canonical Return (CR)

**Independent research lineage, opened 2026-11-22 (America/Los_Angeles) / 2026-11-23 UTC.** This Space is deliberately separate from monocap-tech/weil-lab, monocap-tech/weil, the original RH/provenance corpus, and RH-HODGE. No inherited proof acceptance, workflow trigger, ledger mutation, or mathematical claim promotion. Source repositories are read-only comparators unless separately authorized.

## Motivating question

Can the exact structures and constraints reconstructed by Weil/Weil-Lab—signed arithmetic pairings, phase-sensitive prime-power couplings, finite Gram forms, support restrictions, operator/domain interfaces, and proved obstructions—**recover, constrain, or extend** the canonical constructions that motivated the original September 2026 RH research? A recovered identity is not a new result; a formal dictionary is not an equivalence unless proved; an identified obstruction is a legitimate outcome.

## Author/order

1. **CR-SZ — Masatoshi Suzuki**, screw-function continuous kernels, localized Weil operators, and generalized eigenvalue formulation (active).

2. **CR-MO — Masanori Morishita**, Deninger ↔ Connes–Consani phase-space/orbit map (queued).

3. **CR-CC — Alain Connes / Caterina Consani** (interpreting the user's “Kanes” provisionally as Connes; include CCM where appropriate), Sonin/compression and adelic constructions (queued).

4. **CR-DE — Christopher Deninger**, flow, conjectural arithmetic cohomology, Hodge metric/spectral rigidity (queued).

This order is investigative, not a hierarchy of validity. If “Kanes” denotes a different author, update only that entry after source identification.

## Every investigation must report

- **Primary source/version and exact theorem/equation location**; frozen research source/commit where relevant.

- **Source and target objects**, including common test domain, anti-/sesquilinearity convention, signs, both zeta poles, complete Gamma part, arithmetic weights, scaling and normalizations.

- **Candidate map and equality**, clearly separating full-form equality from termwise reconstruction; operator domains, adjoints, completions and boundary conditions must be checked.

- **Invariants and falsifiers**, including the IV030 phase triple (1, −1, i), support-extension collision graph, correct native metric, and source-admissibility.

- **Outcome class:** RECOVERY (already canonical), OBSTRUCTION, NEW CONSEQUENCE, or OPEN. Record explicitly whether a theorem is only analytic, computationally checked, or Lean-accepted.

- **Potential RH-repository consumer** (STRICT marked-observability, BAL exact-gap rigidity, or another precisely typed obligation), *without* activating a stopped lineage unless a carrier-preserving theorem is proved.

## Baseline and boundaries

- F40 numerical HD/A gate custody: **916/968 accepted**, as of NF781 independent receipt (2026-11-23 UTC). This is not a complete physical matrix or global positivity theorem.

- F41 physical Gram/congruence completion: open. Original unrestricted Weil positivity and RH: open.

- IV030 source-level phase obstruction is an explicit analytic construction with exact rational checks, **not a Lean-accepted new theorem**.

- Original RH September regroup ends in STRICT/BAL theorem deficits; PRED-SEAM and HODGE contain nontransfer determinations. Retain those as preconditions, not as invitations to resume stopped searches automatically.

## First primary references

- Suzuki, *Weil’s quadratic form via the screw function*, [arXiv:2606.09096v4](https://arxiv.org/html/2606.09096), version November 15, 2026. Especially (1.3), (8.7), (8.8).

- Suzuki, *Spectral convergence associated with Weil's quadratic form*, [arXiv:2610.32704](https://arxiv.org/abs/2610.32704), **assuming RH** for spectral convergence.

- Morishita, *On a relation between Deninger's foliated dynamical systems and Connes–Consani's adelic spaces*, [arXiv:2508.15971](https://arxiv.org/abs/2508.15971).

- Weil-Lab [IV030 source-level phase no-go](https://github.com/monocap-tech/weil-lab/blob/1061a331c195e30e85dd8eb35a6eff1c794cf054/docs/IV030_PHASE_IV_ORIGINAL_PRIME_PHASE_FIDELITY_NO_GO.md).

- Original RH [September 13 regroup](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_GLOBAL_REGROUP_20260913.md) and [September 14 RH-HODGE handoff](https://github.com/monocap-tech/provenance/blob/main/rh/handoffs/RH_HODGE_PROJECT_HANDOFF_20260914.txt).

## Status

**CR-SZ-001:** arithmetic prime-kink/relative-phase calibration on Suzuki's continuous kernel — analytic smooth-core recovery; **NOT Lean**.

**CR-SZ-002:** full scalar signed Weil formula recovery **PASS analytically**: both zeta pole cross-moments, signed prime sum, Gamma/digamma normalization and Fourier-scale/argument conventions match exactly on smooth test functions. **Source-carrier/completed-domain equivalence and new Lean kernel acceptance remain OPEN.** [Detailed analysis](https://chatgpt.com/space/page_8c63ddb9ebb481919002432bc72dab3b).

**CR-SZ-003:** **analytic COMPLETE-FORM phase-splitting bound established** on nonnegative translated smooth bump pairs, including original two poles, Gamma, and n=9 prime term. Exact admissibility threshold \(a>\log3\); normalized full mixed coefficient \(-c/J>83/250\), uniform over selected narrow bumps; 2×2 generalized-spectrum gap \(>83/125\); only a \(q\)-dependent **upper bound** for full Suzuki \(\lambda_a\). This extends/recovers related PI027 phase-dominance estimates; no new RH result or Lean compilation. [Full CR-SZ-003 analysis](https://chatgpt.com/space/page_f22d2c1da1f0819181937662f715a667).

**CR-MO-001 (Morishita):** [Theorem 3.6 / orbitwise unitary and metric-transfer obstruction](https://chatgpt.com/space/page_b6e89919038c8191bfbdd6d202544491). Verified the v5 primary paper. Each matched primitive prime circle is carried by time reversal and thus is an orientation-reversing circle homeomorphism. Normalized *orbitwise* Haar pullback is unitary and reverses Fourier mode sign, while no **global** Hilbert/Hodge/Weil pairing is transported. Elementary consequence, not an RH or Lean result; original PRED no-transfer persists.

**CR-CC-001:** [Full window/trace source dictionary, explicit pole-null four-cell phase family and conditional Schur gate](https://chatgpt.com/space/page_9bd218795c188191b89f97a8f592b814). Using Connes–Consani 2020, Connes–Consani–Moscovici 2024 and Andrews–Connes–La Gatta 2026 (Nov 9), the smooth original signed Weil form maps exactly into the multiplicative window under \(u=e^x\). The old IV030 two-bump span has no nonzero two-pole-null member; an explicit four-bump annihilator supplies real and complex pole-null phase probes with identical support/mass and distinct genuine prime terms. CCM finite-place Sonin Hilbert isomorphism is generally nonunitary in the native metric, with explicit finite-place product bounds. The new 2026 conditional Schur reduction still requires (H2'), (R1), (R2) and a uniform transition/pole Schur lower bound; no RH or Lean claim.

**CR-DE-001:** [Deninger Hodge weight, finite Lyapunov test and affine time-reversal obstruction](https://chatgpt.com/space/page_6d0eb966bbd48191b68024af7a4c4bef). Directly inspected Deninger's arXiv:2204.02714v2, properties 2.1–2.7 and Theorem 2.13, plus the original RH PRED-SEAM-16 and HODGE-R133 STOP receipts. **Analytic compatibility theorem:** two positive weight-one Hermitian generators cannot be linked by an *isometric* unshifted flow-reversal map \(\Theta_AT=-T\Theta_D\); a lawful centered-flow reversal must satisfy \(\Theta_AT=T(I-\Theta_D)\). An indefinite-metric 2×2 model demonstrates the weight-one identity alone does not exclude off-line spectral pairs. The Hodge structure remains conjectural, and no original RH or Lean theorem is promoted.

**CR-INTEGRATION-001** — [Five-invariant compatibility audit, native Gram-pencil congruence and explicit finite-prime Hodge-weight model](https://chatgpt.com/space/page_87cd083515c48191b13fe544e63dfc1e). Outcome: the four canonical approaches share **selected** exact formal data but there is **no established global Weil/Hodge/spectral equivalence**. Standard Hermitian-pencil congruence requires transport of both quadratic form and positive native metric. The CCM finite-prime Mellin factor has positive **even** modulus square; on a doubled ambient Mellin Hilbert space it supports an exact algebraic weight-one alternating pairing, positive Hermitian metric, quaternionic \(J^2=-I\), and affine time reversal. This is a **toy formal compatibility model**, NOT arithmetic cohomology or an RH theorem. The selected next noncircular theorem is **CR-INT-T1A**: verify complex polarization and trace-class admissibility of Connes's pole-null semilocal trace on the explicit real four-cell subspace, then attempt a full mixed Gram-pencil attachment. No Lean/new RH result.

**CR-INT-T1A — [Complex mixed semilocal trace and full-Weil pole-null phase ordering](https://chatgpt.com/space/page_b7eef8c3a0dc8191b6e8441b7e0b53c3)**. Reopened Andrews–Connes–La Gatta arXiv:2611.13129v1 §8.2, Theorem 8.1, Lemma 8.3, and its cited full module-Schwartz trace theorem. **Analytically passed** the mixed trace-class domain for complex compact smooth \(f,g\) by a quantized-commutator factorization and established the exact full signed Weil/pole/trace identity on the two-Mellin-pole-null subspace. The prior four-bump polynomial family \(\beta=1,i,-1\) has the *same support, unit mass, and interaction graph*; its **full original Weil form** and semilocal **trace** values are strictly ordered, with each successive gap >1/30, despite all Gamma terms retained. This is a restricted analytic source-attachment consequence of existing canonical trace formulas, **not Lean, not global positivity, not an established novel research theorem**.

**CR-INT-T1B — [Exact three-route pole-null 2×2 Gram pencil and Gamma reinforcement](https://chatgpt.com/space/page_6acd44207bf08191b5e7607534de674d).** Exact real pole-null generators \(e_0=\phi_0-s\phi_1+\phi_2\), \(e_1=\phi_1-s\phi_2+\phi_3\), \(s=4/\sqrt3\), have the same nonorthogonal physical mass \(M=J\bigl[\begin{smallmatrix}22/3&-8\sqrt3/3\\-8\sqrt3/3&22/3\end{smallmatrix}\bigr]\) in the original, Connes semilocal and Suzuki derivative-kernel presentations. The complete signed Weil numerator \(H\) was derived exactly in terms of single-bump Gamma self contribution \(g_0\), cross integrals \(A_{1,2,3}\), and genuine \(3,9,27\) prime powers; an independent SymPy exact audit PASSED. The two selected generalized eigenvalues obey **\(\mu_-−\mu_+>10\sqrt3\log3/73\approx0.260665\)** because the *actual Gamma term reinforces*, rather than opposes, this family’s prime phase splitting (via \(A_2<A_1/\sqrt3\)). The normalized \(1,i,-1\) complete Weil values have minimum adjacent gap greater than \((55\sqrt3-60)\log3/803\approx0.048244\). The pole-null finite-source semilocal trace form is invariant under adding places beyond the window’s prime threshold, **without implying a global Hilbert/operator isomorphism**. This is an analytic finite-compression theorem, not new Lean or RH, and no full-spectrum positivity/sign follows.

**CR-INT-T1C — [Exact rational Gamma enclosures and positive two-dimensional pole-null pencil](https://chatgpt.com/space/page_4092c1d91d848191ad6bc3aecbf7d2d6).** Fixed an actual C∞ nonnegative triangle–mollifier convolution with width 1/10000 and support radius a=2, and the same four-cell/two-Mellin-pole-null family. The executed pure-rational interval auditor yields g0/J ∈ [6.79514722,13.44129808], A1,2,3/J rigorous positive intervals, and absolute selected generalized eigenvalues μ+ ∈ [7.28187642,13.92804353], μ− ∈ [7.54256106,14.18874262]. Thus **the entire selected 2D full signed Weil restriction is strictly positive**, and its correlated eigenvalue splitting lies in [0.26068352,0.26070960]. No quadrature-error guess was used; all decimal endpoints were rationally rounded outward. **This is not a whole-operator lower bound, new Lean certificate, or RH claim.**

**Next: CR-INT-T1D** — relax support geometry and test how the finite pole-null minimum behaves, while separately identifying the missing whole-complement coercivity hypothesis.
