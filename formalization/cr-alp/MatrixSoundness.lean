import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace CRALP

open scoped BigOperators

def quad {n : ℕ} (A : Fin n → Fin n → ℝ) (x : Fin n → ℝ) : ℝ :=
  ∑ i, ∑ j, A i j * (x i * x j)

theorem entryLower (e x y ε : ℝ) (he : |e| ≤ ε) (hε : 0 ≤ ε) :
    -(ε / 2) * (x ^ 2 + y ^ 2) ≤ e * (x * y) := by
  have hxy : |x * y| ≤ (x ^ 2 + y ^ 2) / 2 := by
    rw [abs_le]
    constructor <;> nlinarith [sq_nonneg (x+y), sq_nonneg (x-y)]
  have h : |e * (x * y)| ≤ ε * ((x ^ 2 + y ^ 2) / 2) := by
    rw [abs_mul]
    exact mul_le_mul he hxy (abs_nonneg _) hε
  have hneg := neg_abs_le (e * (x * y))
  calc
    -(ε / 2) * (x ^ 2 + y ^ 2) =
        -(ε * ((x ^ 2 + y ^ 2) / 2)) := by ring
    _ ≤ -|e * (x * y)| := neg_le_neg h
    _ ≤ e * (x * y) := hneg

theorem uniformErrorLower {n : ℕ} (R : Fin n → Fin n → ℝ)
    (x : Fin n → ℝ) (ε : ℝ) (hε : 0 ≤ ε)
    (hR : ∀ i j, |R i j| ≤ ε) :
    -(ε * (n : ℝ)) * (∑ i, x i ^ 2) ≤ quad R x := by
  unfold quad
  calc
    -(ε * (n : ℝ)) * (∑ i, x i ^ 2) =
        ∑ i : Fin n, ∑ j : Fin n, -(ε / 2) * (x i ^ 2 + x j ^ 2) := by
      simp only [mul_add, Finset.sum_add_distrib, Finset.sum_const,
        Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, ← Finset.mul_sum]
      ring
    _ ≤ ∑ i : Fin n, ∑ j : Fin n, R i j * (x i * x j) := by
      exact Finset.sum_le_sum (fun i _ =>
        Finset.sum_le_sum (fun j _ => entryLower _ _ _ _ (hR i j) hε))

theorem gramQuadratic {n : ℕ} (L : Fin n → Fin n → ℝ)
    (x : Fin n → ℝ) :
    quad (fun i j => ∑ k, L i k * L j k) x =
      ∑ k, (∑ i, L i k * x i) ^ 2 := by
  unfold quad
  calc
    (∑ i : Fin n, ∑ j : Fin n, (∑ k : Fin n, L i k * L j k) * (x i * x j)) =
        ∑ i : Fin n, ∑ j : Fin n, ∑ k : Fin n, (L i k * x i) * (L j k * x j) := by
      apply Finset.sum_congr rfl
      intro i _
      apply Finset.sum_congr rfl
      intro j _
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro k _
      ring
    _ = ∑ i : Fin n, ∑ k : Fin n, ∑ j : Fin n, (L i k * x i) * (L j k * x j) := by
      apply Finset.sum_congr rfl
      intro i _
      exact Finset.sum_comm
    _ = ∑ k : Fin n, ∑ i : Fin n, ∑ j : Fin n, (L i k * x i) * (L j k * x j) :=
      Finset.sum_comm
    _ = ∑ k : Fin n, (∑ i : Fin n, L i k * x i) ^ 2 := by
      apply Finset.sum_congr rfl
      intro k _
      rw [pow_two, Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro i _
      rw [Finset.mul_sum]

theorem diagonalQuadratic {n : ℕ} (δ : ℝ) (x : Fin n → ℝ) :
    quad (fun i j => if i = j then δ else 0) x = δ * (∑ i, x i ^ 2) := by
  unfold quad
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  simp [ite_mul, pow_two]

/-- A real Gram-plus-small-residual witness implies a strict lower bound.
    The entrywise witness premises are separate exact arithmetic obligations. -/
theorem gramResidualLower {n : ℕ} (A L R : Fin n → Fin n → ℝ)
    (δ ε : ℝ) (hε : 0 ≤ ε)
    (hR : ∀ i j, |R i j| ≤ ε)
    (hA : ∀ i j, A i j = (∑ k, L i k * L j k) +
      (if i = j then δ else 0) + R i j)
    (x : Fin n → ℝ) :
    (δ - ε * (n : ℝ)) * (∑ i, x i ^ 2) ≤ quad A x := by
  have hsplit : quad A x =
      quad (fun i j => ∑ k, L i k * L j k) x +
      quad (fun i j => if i = j then δ else 0) x + quad R x := by
    simp only [quad, hA, add_mul, Finset.sum_add_distrib]
  rw [hsplit, gramQuadratic, diagonalQuadratic]
  have hg : 0 ≤ ∑ k : Fin n, (∑ i : Fin n, L i k * x i) ^ 2 :=
    Finset.sum_nonneg (fun _ _ => sq_nonneg _)
  have hr := uniformErrorLower R x ε hε hR
  rw [sub_mul]
  linarith

end CRALP

#print axioms CRALP.entryLower
#print axioms CRALP.uniformErrorLower
#print axioms CRALP.gramQuadratic
#print axioms CRALP.diagonalQuadratic
#print axioms CRALP.gramResidualLower
