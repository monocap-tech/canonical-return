# CR-ALP LEAN-02A-R1 — isolated double-sum rewrite failure

Incident CR-ALP/LEAN-02AR1/20261127/06; update to incident 05, distinct diagnostic. Run 46088611748, job 140803965202, activation 07f4deab7d23069ceb65f7624acb8cc94a6ff2c1, source blob 5d8aad0bd50c68ec9d9ae4109d94d168dc6008ea. All six external exact matrix witnesses passed, setup passed. Native compilation failed at line 41: reverse rewrite of Finset.mul_sum does not match the nested expression. The hsum equality itself elaborated. entryLower, gramQuadratic, diagonalQuadratic printed only permitted axioms; uniformErrorLower and gramResidualLower depend on sorryAx due to the failed rewrite. No accepted native module or durable result push. The first repair fixed the missing import; new diagnostic is not an unchanged retry.

Both failed archives independently authenticated locally: original 15240146393 SHA256 daf7a85c6b3503a4a8b93441551a311ca905e1876699d7d77e22321796b5c9b7; R1 15240890776 SHA256 3b1051cd4ec4f9d935446bc2f9d212037f0df6683ec704be2fdc595261ff2636. ZIP CRC PASS both. Both witnesses.json members byte-identical (1003035 bytes, SHA256 d729ae089f673a2426cfe0e29934cf88ec2074f30fbbd8ca0614f37bd236ee57). R1 axioms.log 1104 bytes SHA256 69963efbba2659c5eeeae6cdcd19f4782f59c2bbb62abda0e2c5a3d49d4fabec. Full raw archives remain local and expiring in Actions, not permanently committed. Central indexing pending. No matrix witness or inequality failed.

Bounded repair R2: replace the fragile reverse rewrite with an explicit intermediate equality to the scalar times the unweighted double sum. First use hsum and ring; then distribute the scalar into both sums using forward Finset.mul_sum via simp only. All five statements and strict verification criteria unchanged. Retain original source and failure records. If the same failure recurs, stop that proof path for reassessment rather than retry blindly. No retirement credit yet.

```text
2026-11-28T00:41:07.7346088Z MatrixSoundnessR1.lean:41:10: error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
2026-11-28T00:41:07.7346813Z   ∑ i ∈ ?s, ?a * ?f i
2026-11-28T00:41:07.7347051Z in the target expression
2026-11-28T00:41:07.7347640Z   -(ε * ↑n) * ∑ i, x i ^ 2 = ∑ i, ∑ j, -(ε / 2) * (x i ^ 2 + x j ^ 2)
2026-11-28T00:41:07.7347909Z 
2026-11-28T00:41:07.7348011Z n : ℕ
2026-11-28T00:41:07.7348225Z R : Fin n → Fin n → ℝ
2026-11-28T00:41:07.7348464Z x : Fin n → ℝ
2026-11-28T00:41:07.7348670Z ε : ℝ
2026-11-28T00:41:07.7348910Z hε : 0 ≤ ε
2026-11-28T00:41:07.7349153Z hR : ∀ (i j : Fin n), |R i j| ≤ ε
2026-11-28T00:41:07.7349882Z hsum : ∑ i, ∑ j, (x i ^ 2 + x j ^ 2) = ↑n * ∑ i, x i ^ 2 + ↑n * ∑ i, x i ^ 2
2026-11-28T00:41:07.7350421Z ⊢ -(ε * ↑n) * ∑ i, x i ^ 2 = ∑ i, ∑ j, -(ε / 2) * (x i ^ 2 + x j ^ 2)
2026-11-28T00:41:07.7350947Z 'CRALP.entryLower' depends on axioms: [propext, Classical.choice, Quot.sound]
2026-11-28T00:41:07.7351613Z 'CRALP.uniformErrorLower' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:41:07.7352509Z 'CRALP.gramQuadratic' depends on axioms: [propext, Classical.choice, Quot.sound]
2026-11-28T00:41:07.7353163Z 'CRALP.diagonalQuadratic' depends on axioms: [propext, Classical.choice, Quot.sound]
2026-11-28T00:41:07.7354096Z 'CRALP.gramResidualLower' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:41:07.8419455Z ##[error]Process completed with exit code 1.
2026-11-28T00:41:07.8511425Z ##[group]Run actions/upload-artifact@v4
2026-11-28T00:41:07.8511736Z with:
2026-11-28T00:41:07.8511956Z   name: CR-ALP-LEAN-02A-R1-matrix-soundness
2026-11-28T00:41:07.8512292Z   path: /tmp/cr-alp-lean02ar1-evidence/
2026-11-28T00:41:07.8512585Z   if-no-files-found: warn
2026-11-28T00:41:07.8512825Z   compression-level: 6
2026-11-28T00:41:07.8513045Z   overwrite: false
```
