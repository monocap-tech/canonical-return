import RationalTrace

namespace CRALP.PrimeLogTrace
open Finset
open scoped BigOperators
open CRALP.DirectedInterval CRALP.DirectedEvaluation CRALP.RationalTrace

def termNum (p i : ℕ) : ℤ := (p^(2*i+1) : ℕ)
def termDen (q i : ℕ) : ℤ := (q^(2*i+1)*(2*i+1) : ℕ)

theorem termIdentity (p q i : ℕ) :
    (termNum p i : ℝ)/(termDen q i : ℝ) =
      ((p:ℝ)/q)^(2*i+1)/(2*i+1 : ℕ) := by
  simp only [termNum,termDen,Int.cast_mul,Int.cast_pow,Int.cast_natCast,
    Nat.cast_mul,Nat.cast_pow,div_pow,div_div]

theorem reducedTrace (x : ℝ) (p q : ℕ) (lo hi : ℕ → ℤ) (L U : ℤ)
    (hx : 0 < x) (hm : (x-1)/(x+1) = (p:ℝ)/q)
    (hr : |(x-1)/(x+1)| ≤ 7/20)
    (hrow : ∀ i : Fin 400, row (termNum p) (termDen q) lo hi i)
    (hl : L=2*(∑ i ∈ range 400, lo i)-1)
    (hu : U=2*(∑ i ∈ range 400, hi i)+1) :
    Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.log x) := by
  have h := sumTrace (termNum p) (termDen q) lo hi hrow
  simp only [termIdentity] at h
  have hp : Contains
      (2*((∑ i ∈ range 400, (lo i:ℝ))/logScale))
      (2*((∑ i ∈ range 400, (hi i:ℝ))/logScale))
      (ElementarySeries.logPoly ((x-1)/(x+1))) := by
    rw [hm]
    unfold ElementarySeries.logPoly
    exact ⟨mul_le_mul_of_nonneg_left h.1 (by norm_num),
      mul_le_mul_of_nonneg_left h.2 (by norm_num)⟩
  have ha := logAttach hx hr hp
  have hL : (L:ℝ)/logScale =
      2*((∑ i ∈ range 400, (lo i:ℝ))/logScale)-ElementarySeries.unit := by
    rw [hl,unitEq]
    push_cast <;> ring
  have hU : (U:ℝ)/logScale =
      2*((∑ i ∈ range 400, (hi i:ℝ))/logScale)+ElementarySeries.unit := by
    rw [hu,unitEq]
    push_cast <;> ring
  rw [hL,hU]
  exact ha

/-- Nonnegative dyadic reductions used by the literal prime-window logs. -/
theorem scaleTrace (x y : ℝ) (k : ℕ) (l u l2 u2 L U : ℤ)
    (hx : 0 < x) (hy : x/(2:ℝ)^k = y)
    (hr : Contains ((l:ℝ)/logScale) ((u:ℝ)/logScale) (Real.log y))
    (h2 : Contains ((l2:ℝ)/logScale) ((u2:ℝ)/logScale) (Real.log 2))
    (hl : L=l+(k:ℤ)*l2) (hu : U=u+(k:ℤ)*u2) :
    Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.log x) := by
  have he : Real.log x = Real.log y+(k:ℝ)*Real.log 2 := by
    rw [← hy,Real.log_div hx.ne' (by positivity : (2:ℝ)^k ≠ 0),Real.log_pow]
    ring
  have hk : (0:ℝ) ≤ k := by positivity
  have hs : Contains ((k:ℝ)*((l2:ℝ)/logScale))
      ((k:ℝ)*((u2:ℝ)/logScale)) ((k:ℝ)*Real.log 2) :=
    ⟨mul_le_mul_of_nonneg_left h2.1 hk,mul_le_mul_of_nonneg_left h2.2 hk⟩
  have h := DirectedInterval.add hr hs
  have hL : (L:ℝ)/logScale =
      (l:ℝ)/logScale+(k:ℝ)*((l2:ℝ)/logScale) := by
    rw [hl]
    push_cast <;> ring
  have hU : (U:ℝ)/logScale =
      (u:ℝ)/logScale+(k:ℝ)*((u2:ℝ)/logScale) := by
    rw [hu]
    push_cast <;> ring
  rw [hL,hU,he]
  exact h

end CRALP.PrimeLogTrace
#print axioms CRALP.PrimeLogTrace.termIdentity
#print axioms CRALP.PrimeLogTrace.reducedTrace
#print axioms CRALP.PrimeLogTrace.scaleTrace
