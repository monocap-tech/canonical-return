import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

/-! First native CR-ALP retirement component: scalar inequalities.
The physical coefficient bounds and the infinite-dimensional receiver
attachment are explicit hypotheses, not imported proof declarations. -/

namespace CRALP

theorem directSchurLower (hmin corrected defect schur rho : ℝ)
    (hc : hmin ≤ corrected) (hd : defect ≤ rho * hmin)
    (heq : schur = corrected - defect) :
    (1 - rho) * hmin ≤ schur := by
  nlinarith

theorem squareCeiling {x a : ℝ} (hx : 0 ≤ x) (ha : 0 ≤ a)
    (h : x ^ 2 ≤ a ^ 2) : x ≤ a := by
  by_contra hn
  have hlt : a < x := lt_of_not_ge hn
  have hp : 0 < (x - a) * (x + a) :=
    mul_pos (sub_pos.mpr hlt) (by linarith)
  nlinarith

/-- The normalized inverse-response loss, with rational ceilings on both
    square roots. The hypothesis hphysical is the remaining physical receiver
    contract, not a conclusion of the finite arithmetic. -/
theorem poleBudget {E sigma u x y uMax xMax yMax loss : ℝ}
    (hE : 0 < E) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hxMax : 0 ≤ xMax) (hyMax : 0 ≤ yMax)
    (hu : u ≤ uMax) (hxs : x ^ 2 ≤ xMax ^ 2)
    (hys : y ^ 2 ≤ yMax ^ 2)
    (hbudget : uMax + (xMax + yMax) ^ 2 ≤ loss)
    (hphysical : E * (1 - u - (x + y) ^ 2) ≤ sigma) :
    E * (1 - loss) ≤ sigma := by
  have hxl := squareCeiling hx hxMax hxs
  have hyl := squareCeiling hy hyMax hys
  have hsum : x + y ≤ xMax + yMax := add_le_add hxl hyl
  have hp : 0 ≤ (xMax + yMax - (x + y)) *
      (xMax + yMax + (x + y)) :=
    mul_nonneg (sub_nonneg.mpr hsum) (by linarith)
  have hs : (x + y) ^ 2 ≤ (xMax + yMax) ^ 2 := by
    nlinarith
  have hl : u + (x + y) ^ 2 ≤ loss := by linarith
  exact le_trans
    (mul_le_mul_of_nonneg_left (show 1 - loss ≤ 1 - u - (x + y) ^ 2 by linarith) hE.le)
    hphysical

/-- AQ2 rounded normalized budget, cosh. -/
theorem coshTerminal {E sigma u x y : ℝ}
    (hE : 0 < E) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hu : u ≤ 425975815994492321 / 1000000000000000000)
    (hb : x ^ 2 ≤ 277101 / 250000000000000000000000000000)
    (hz : y ^ 2 ≤ (251443 / 248557) * u)
    (hphysical : E * (1 - u - (x + y) ^ 2) ≤ sigma) :
    E * (1 - 856897646 / 1000000000) ≤ sigma := by
  have hxs : x ^ 2 ≤ (131601 / 125000000000000000 : ℝ) ^ 2 :=
    le_trans hb (by norm_num)
  have hys : y ^ 2 ≤ (656446364299909703 / 1000000000000000000 : ℝ) ^ 2 :=
    le_trans hz (le_trans (mul_le_mul_of_nonneg_left hu (by norm_num : (0 : ℝ) ≤ 251443 / 248557)) (by norm_num))
  exact poleBudget hE hx hy (by norm_num) (by norm_num) hu hxs hys (by norm_num) hphysical

/-- AQ2 rounded normalized budget, sinh. -/
theorem sinhTerminal {E sigma u x y : ℝ}
    (hE : 0 < E) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hu : u ≤ 188882199696291911 / 500000000000000000)
    (hb : x ^ 2 ≤ 170643 / 1000000000000000000000000000000)
    (hz : y ^ 2 ≤ (108431 / 141569) * u)
    (hphysical : E * (1 - u - (x + y) ^ 2) ≤ sigma) :
    E * (1 - 667102967 / 1000000000) ≤ sigma := by
  have hxs : x ^ 2 ≤ (41309 / 100000000000000000 : ℝ) ^ 2 :=
    le_trans hb (by norm_num)
  have hys : y ^ 2 ≤ (107580401001730913 / 200000000000000000 : ℝ) ^ 2 :=
    le_trans hz (le_trans (mul_le_mul_of_nonneg_left hu (by norm_num : (0 : ℝ) ≤ 108431 / 141569)) (by norm_num))
  exact poleBudget hE hx hy (by norm_num) (by norm_num) hu hxs hys (by norm_num) hphysical

theorem bothMarginsPositive (Ecosh Esinh sc ss : ℝ)
    (he : 0 < Ecosh) (ho : 0 < Esinh)
    (hc : Ecosh * (1 - 856897646 / 1000000000) ≤ sc)
    (hs : Esinh * (1 - 667102967 / 1000000000) ≤ ss) :
    0 < sc ∧ 0 < ss := by
  exact ⟨lt_of_lt_of_le (mul_pos he (by norm_num)) hc,
    lt_of_lt_of_le (mul_pos ho (by norm_num)) hs⟩

end CRALP

#print axioms CRALP.directSchurLower
#print axioms CRALP.squareCeiling
#print axioms CRALP.poleBudget
#print axioms CRALP.coshTerminal
#print axioms CRALP.sinhTerminal
#print axioms CRALP.bothMarginsPositive
