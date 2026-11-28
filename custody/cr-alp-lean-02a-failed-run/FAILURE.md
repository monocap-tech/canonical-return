# CR-ALP LEAN-02A — failed native matrix soundness verification

Incident CR-ALP/LEAN-02A/20261127/05 (UTC 2026-11-28). Run 46087901668, job 140801706699, activation 5597e97fde574fd810a14a4200a6e0a64d75c952. Source MatrixSoundness.lean blob dd6dcfe96f53be5c3b5d8138d910a957b96a6ab6. Pinned Lean 4.34.0 and Mathlib 5ed2965256430c3649e86755f9576b54eca72435. All six exact external Gram residual witness tests passed and isolated setup succeeded. Native compile failed at missing Finset.mul_sum / Finset.sum_mul and resulting unresolved sum goals; entryLower alone printed an allowed axiom set. The other four declarations printed sorryAx due to elaboration errors, not valid proofs. No accepted native result, no durable runner result push. No mathematical matrix positivity failure.

Artifact 15240146393, 305486 bytes, API SHA256 daf7a85c6b3503a4a8b93441551a311ca905e1876699d7d77e22321796b5c9b7. Downloaded archive file CR-ALP-LEAN-02A-failed.zip, local only; byte authentication remains separate. Raw artifact still expiring; this record preserves exact diagnostic excerpt and identity, not entire raw archive. Central incident indexing pending.

At the exact pinned Mathlib revision, Mathlib/Algebra/BigOperators/Ring/Finset.lean (blob bedfc7ff7ff4932e5fd1979ccb663071cc3664fc) contains both missing lemmas. Bounded repair imports that module and separates the double-sum factorization from scalar normalization. All proof statements, matrix inputs, witness acceptance conditions and strict named axiom checks retained. Prevention: verify exact lemma imports and compile generic sum lemmas before data-heavy execution. No matrix data kernel check or retirement promotion until successful native evidence.

```text
2026-11-28T00:31:17.9868153Z MatrixSoundness.lean:2:0: warning: 
2026-11-28T00:31:17.9902027Z 'Mathlib.Data.Real.Basic' has been deprecated: please replace this import by
2026-11-28T00:31:17.9931535Z 
2026-11-28T00:31:17.9941896Z import Mathlib.Basic.Real.Basic
2026-11-28T00:31:17.9942455Z MatrixSoundness.lean:37:60: error(lean.unknownIdentifier): Unknown constant `Finset.mul_sum`
2026-11-28T00:31:17.9943033Z Try this:
2026-11-28T00:31:17.9943314Z   [apply] ring_nf
2026-11-28T00:31:17.9943596Z   
2026-11-28T00:31:17.9944009Z   The `ring` tactic failed to close the goal. Use `ring_nf` to obtain a normal form.
2026-11-28T00:31:17.9944469Z     
2026-11-28T00:31:17.9945764Z   Note that `ring` works primarily in *commutative* rings. If you have a noncommutative ring, abelian group or module, consider using `noncomm_ring`, `abel` or `module` instead.
2026-11-28T00:31:17.9946453Z MatrixSoundness.lean:35:68: error: unsolved goals
2026-11-28T00:31:17.9947017Z n : ℕ
2026-11-28T00:31:17.9947322Z R : Fin n → Fin n → ℝ
2026-11-28T00:31:17.9947661Z x : Fin n → ℝ
2026-11-28T00:31:17.9947951Z ε : ℝ
2026-11-28T00:31:17.9948483Z hε : 0 ≤ ε
2026-11-28T00:31:17.9948799Z hR : ∀ (i j : Fin n), |R i j| ≤ ε
2026-11-28T00:31:17.9949407Z ⊢ -(ε * ↑n * ∑ x_1, x x_1 ^ 2) = ↑n * ∑ x_1, ε * x x_1 ^ 2 * (-1 / 2) + ∑ x_1, ε * ↑n * x x_1 ^ 2 * (-1 / 2)
2026-11-28T00:31:17.9950315Z MatrixSoundness.lean:55:10: error(lean.unknownIdentifier): Unknown constant `Finset.sum_mul`
2026-11-28T00:31:17.9950902Z MatrixSoundness.lean:50:80: error: unsolved goals
2026-11-28T00:31:17.9951492Z n : ℕ
2026-11-28T00:31:17.9951849Z L : Fin n → Fin n → ℝ
2026-11-28T00:31:17.9952289Z x : Fin n → ℝ
2026-11-28T00:31:17.9952561Z i : Fin n
2026-11-28T00:31:17.9952874Z a✝¹ : i ∈ Finset.univ
2026-11-28T00:31:17.9953106Z j : Fin n
2026-11-28T00:31:17.9953404Z a✝ : j ∈ Finset.univ
2026-11-28T00:31:17.9953772Z ⊢ (∑ k, L i k * L j k) * (x i * x j) = ∑ k, L i k * x i * (L j k * x j)
2026-11-28T00:31:17.9954293Z MatrixSoundness.lean:68:19: error(lean.unknownIdentifier): Unknown constant `Finset.sum_mul`
2026-11-28T00:31:17.9954857Z MatrixSoundness.lean:65:55: error: unsolved goals
2026-11-28T00:31:17.9955216Z n : ℕ
2026-11-28T00:31:17.9955531Z L : Fin n → Fin n → ℝ
2026-11-28T00:31:17.9955817Z x : Fin n → ℝ
2026-11-28T00:31:17.9956018Z k : Fin n
2026-11-28T00:31:17.9956273Z a✝ : k ∈ Finset.univ
2026-11-28T00:31:17.9956610Z ⊢ ∑ i, ∑ j, L i k * x i * (L j k * x j) = (∑ i, L i k * x i) * ∑ i, L i k * x i
2026-11-28T00:31:17.9957040Z MatrixSoundness.lean:76:6: error(lean.unknownIdentifier): Unknown constant `Finset.mul_sum`
2026-11-28T00:31:17.9957454Z MatrixSoundness.lean:74:71: error: unsolved goals
2026-11-28T00:31:17.9957727Z n : ℕ
2026-11-28T00:31:17.9957972Z δ : ℝ
2026-11-28T00:31:17.9958180Z x : Fin n → ℝ
2026-11-28T00:31:17.9958536Z ⊢ ∑ i, ∑ j, (fun i j => if i = j then δ else 0) i j * (x i * x j) = δ * ∑ i, x i ^ 2
2026-11-28T00:31:17.9958955Z 'CRALP.entryLower' depends on axioms: [propext, Classical.choice, Quot.sound]
2026-11-28T00:31:17.9959413Z 'CRALP.uniformErrorLower' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:31:17.9959943Z 'CRALP.gramQuadratic' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:31:17.9960429Z 'CRALP.diagonalQuadratic' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:31:17.9960923Z 'CRALP.gramResidualLower' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
2026-11-28T00:31:18.0503708Z ##[error]Process completed with exit code 1.
2026-11-28T00:31:18.0723531Z ##[group]Run actions/upload-artifact@v4
```
