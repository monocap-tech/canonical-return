import MatrixSoundnessR2
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic

namespace CRALP

open scoped BigOperators

/-- Computable signed residual. List.ofFn avoids an opaque finite-set fold
    in concrete kernel reduction. -/
def intResidual {n : ℕ} (M N : Fin n → Fin n → ℤ)
    (D S d p : ℤ) (i j : Fin n) : ℤ :=
  d * S^2 * M i j - d * D * (List.ofFn (fun k => N i k * N j k)).sum -
    (if i = j then p * D * S^2 else 0)

abbrev rowOK {n : ℕ} (M N : Fin n → Fin n → ℤ)
    (D S d p : ℤ) (i : Fin n) : Prop :=
  ∀ j : Fin n,
    -(p * D * S^2) ≤ 2 * (n : ℤ) * intResidual M N D S d p i j ∧
    2 * (n : ℤ) * intResidual M N D S d p i j ≤ p * D * S^2

theorem intWitnessLower {n : ℕ} (M N : Fin n → Fin n → ℤ)
    (D S d p : ℤ) (hn : 0 < n) (hD : 0 < D) (hS : 0 < S)
    (hd : 0 < d) (hp : 0 < p)
    (hcheck : ∀ i, rowOK M N D S d p i)
    (x : Fin n → ℝ) :
    ((p : ℝ) / (2 * d)) * (∑ i, x i ^ 2) ≤
      quad (fun i j => (M i j : ℝ) / D) x := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hD' : (0 : ℝ) < D := by exact_mod_cast hD
  have hS' : (0 : ℝ) < S := by exact_mod_cast hS
  have hd' : (0 : ℝ) < d := by exact_mod_cast hd
  have hp' : (0 : ℝ) < p := by exact_mod_cast hp
  let den : ℝ := (d : ℝ) * D * (S : ℝ)^2
  have hden : 0 < den := by dsimp [den]; positivity
  let δ : ℝ := (p : ℝ) / d
  let ε : ℝ := δ / (2 * n)
  have hδ : 0 < δ := div_pos hp' hd'
  have hε : 0 ≤ ε := div_nonneg hδ.le (by positivity)
  let L : Fin n → Fin n → ℝ := fun i j => (N i j : ℝ) / S
  let R : Fin n → Fin n → ℝ := fun i j =>
    (intResidual M N D S d p i j : ℝ) / den
  have hscale : ε * (2 * (n : ℝ) * den) = (p : ℝ) * D * (S : ℝ)^2 := by
    dsimp [ε, δ, den]
    field_simp [ne_of_gt hd', ne_of_gt hn'] <;> ring
  have hR : ∀ i j, |R i j| ≤ ε := by
    intro i j
    have hl := (hcheck i j).1
    have hu := (hcheck i j).2
    have hl' : -((p : ℝ) * D * (S : ℝ)^2) ≤
        2 * (n : ℝ) * (intResidual M N D S d p i j : ℝ) := by exact_mod_cast hl
    have hu' : 2 * (n : ℝ) * (intResidual M N D S d p i j : ℝ) ≤
        (p : ℝ) * D * (S : ℝ)^2 := by exact_mod_cast hu
    have hfac : 0 < 2 * (n : ℝ) := by positivity
    have htop : (intResidual M N D S d p i j : ℝ) ≤ ε * den := by
      apply le_of_mul_le_mul_left _ hfac
      calc
        2 * (n : ℝ) * (intResidual M N D S d p i j : ℝ) ≤
            (p : ℝ) * D * (S : ℝ)^2 := hu'
        _ = 2 * (n : ℝ) * (ε * den) := by rw [← hscale]; ring
    have hbot : -ε * den ≤ (intResidual M N D S d p i j : ℝ) := by
      apply le_of_mul_le_mul_left _ hfac
      calc
        2 * (n : ℝ) * (-ε * den) = -((p : ℝ) * D * (S : ℝ)^2) := by
          rw [← hscale]; ring
        _ ≤ 2 * (n : ℝ) * (intResidual M N D S d p i j : ℝ) := hl'
    change |(intResidual M N D S d p i j : ℝ) / den| ≤ ε
    rw [abs_le]
    exact ⟨(le_div_iff₀ hden).2 hbot, (div_le_iff₀ hden).2 htop⟩
  have hA : ∀ i j, (M i j : ℝ) / D =
      (∑ k, L i k * L j k) + (if i = j then δ else 0) + R i j := by
    intro i j
    have hsum : (∑ k, L i k * L j k) =
        (∑ k, (N i k : ℝ) * (N j k : ℝ)) / (S : ℝ)^2 := by
      calc
        (∑ k, L i k * L j k) =
            ∑ k, ((N i k : ℝ) * (N j k : ℝ)) / (S : ℝ)^2 := by
          apply Finset.sum_congr rfl
          intro k _
          dsimp [L]
          rw [div_mul_div_comm, pow_two]
        _ = _ := by simp only [div_eq_mul_inv, Finset.sum_mul]
    rw [hsum]
    dsimp [R, intResidual, den, δ]
    rw [Fin.sum_ofFn]
    split_ifs <;> push_cast <;>
      field_simp [ne_of_gt hD', ne_of_gt hS', ne_of_gt hd'] <;> ring
  have hbound := gramResidualLower (fun i j => (M i j : ℝ) / D) L R δ ε hε hR hA x
  have hmargin : δ - ε * (n : ℝ) = (p : ℝ) / (2 * d) := by
    dsimp [ε, δ]
    field_simp [ne_of_gt hd', ne_of_gt hn'] <;> ring
  rwa [hmargin] at hbound

end CRALP

#print axioms CRALP.intWitnessLower
