import DirectedEvaluation

namespace CRALP.RationalTrace
open Finset
open scoped BigOperators
open CRALP.DirectedInterval CRALP.DirectedEvaluation

def row (num den lo hi : ℕ → ℤ) (i : Fin 400) : Prop :=
  0 < den i.val ∧ lo i.val * den i.val ≤ num i.val * logScale ∧
    num i.val * logScale ≤ hi i.val * den i.val

instance (num den lo hi : ℕ → ℤ) (i : Fin 400) : Decidable (row num den lo hi i) := by
  unfold row
  infer_instance

theorem sumTrace (num den lo hi : ℕ → ℤ)
    (hrow : ∀ i : Fin 400, row num den lo hi i) :
    Contains ((∑ i ∈ range 400, (lo i : ℝ))/logScale)
      ((∑ i ∈ range 400, (hi i : ℝ))/logScale)
      (∑ i ∈ range 400, (num i : ℝ)/(den i : ℝ)) := by
  have hs : 0 < logScale := by unfold logScale; decide
  have ht (i : ℕ) (hi' : i ∈ range 400) :
      Contains ((lo i : ℝ)/logScale) ((hi i : ℝ)/logScale)
        ((num i : ℝ)/(den i : ℝ)) := by
    let j : Fin 400 := ⟨i,mem_range.mp hi'⟩
    have h := hrow j
    change 0 < den i ∧ lo i * den i ≤ num i * logScale ∧
      num i * logScale ≤ hi i * den i at h
    exact scaledRatio (lo i) (hi i) (num i) (den i) logScale h.1 hs h.2.1 h.2.2
  constructor
  · rw [Finset.sum_div]
    exact Finset.sum_le_sum (fun i hi' => (ht i hi').1)
  · rw [Finset.sum_div]
    exact Finset.sum_le_sum (fun i hi' => (ht i hi').2)

def atanNum (p i : ℕ) : ℤ := (-1)^i * ((p^(2*i+1) : ℕ) : ℤ)
def atanDen (q i : ℕ) : ℤ := (q^(2*i+1)*(2*i+1) : ℕ)
def expNum (p i : ℕ) : ℤ := (p^i : ℕ)
def expDen (q i : ℕ) : ℤ := (q^i * Nat.factorial i : ℕ)

theorem atanTerm (p q i : ℕ) :
    (atanNum p i : ℝ)/(atanDen q i : ℝ) =
      (-1:ℝ)^i*((p:ℝ)/q)^(2*i+1)/(2*i+1 : ℕ) := by
  simp only [atanNum,atanDen,Int.cast_mul,Int.cast_pow,Int.cast_neg,Int.cast_one,
    Int.cast_natCast,Nat.cast_mul,Nat.cast_pow,div_pow,mul_div_assoc,div_div] <;> ring

theorem expTerm (p q i : ℕ) :
    (expNum p i : ℝ)/(expDen q i : ℝ) =
      ((p:ℝ)/q)^i/(Nat.factorial i : ℝ) := by
  simp only [expNum,expDen,Int.cast_natCast,Nat.cast_mul,Nat.cast_pow,div_pow,div_div]

theorem unitEq : ElementarySeries.unit = 1/(logScale:ℝ) := by
  simp only [ElementarySeries.unit,logScale,Nat.cast_pow,Nat.cast_ofNat,
    Int.cast_pow,Int.cast_ofNat]

theorem widenTrace {l u L U : ℤ} {x : ℝ}
    (h : Contains ((l:ℝ)/logScale-ElementarySeries.unit)
      ((u:ℝ)/logScale+ElementarySeries.unit) x)
    (hl : L=l-1) (hu : U=u+1) :
    Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) x := by
  rw [hl,hu]
  push_cast
  rw [sub_div,add_div,← unitEq]
  exact h

theorem atanTrace (p q : ℕ) (lo hi : ℕ → ℤ) (L U : ℤ)
    (hr : |(p:ℝ)/q| ≤ 1/3)
    (hrow : ∀ i : Fin 400, row (atanNum p) (atanDen q) lo hi i)
    (hl : L=(∑ i ∈ range 400, lo i)-1)
    (hu : U=(∑ i ∈ range 400, hi i)+1) :
    Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.arctan ((p:ℝ)/q)) := by
  have h := sumTrace (atanNum p) (atanDen q) lo hi hrow
  simp only [atanTerm] at h
  have hp : Contains ((∑ i ∈ range 400, (lo i:ℝ))/logScale)
      ((∑ i ∈ range 400, (hi i:ℝ))/logScale)
      (ElementarySeries.atanPoly 400 ((p:ℝ)/q)) := h
  have ha := atanAttach hr hp
  apply widenTrace (l := ∑ i ∈ range 400, lo i) (u := ∑ i ∈ range 400, hi i) _ hl hu
  push_cast
  exact ha

theorem expTrace (p q : ℕ) (lo hi : ℕ → ℤ) (L U : ℤ)
    (hr : |(p:ℝ)/q| ≤ 3)
    (hrow : ∀ i : Fin 400, row (expNum p) (expDen q) lo hi i)
    (hl : L=(∑ i ∈ range 400, lo i)-1)
    (hu : U=(∑ i ∈ range 400, hi i)+1) :
    Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.exp ((p:ℝ)/q)) := by
  have h := sumTrace (expNum p) (expDen q) lo hi hrow
  simp only [expTerm] at h
  have hp : Contains ((∑ i ∈ range 400, (lo i:ℝ))/logScale)
      ((∑ i ∈ range 400, (hi i:ℝ))/logScale)
      (ElementarySeries.expPoly ((p:ℝ)/q)) := h
  have ha := expAttach hr hp
  apply widenTrace (l := ∑ i ∈ range 400, lo i) (u := ∑ i ∈ range 400, hi i) _ hl hu
  push_cast
  exact ha

end CRALP.RationalTrace
#print axioms CRALP.RationalTrace.sumTrace
#print axioms CRALP.RationalTrace.atanTerm
#print axioms CRALP.RationalTrace.expTerm
#print axioms CRALP.RationalTrace.unitEq
#print axioms CRALP.RationalTrace.widenTrace
#print axioms CRALP.RationalTrace.atanTrace
#print axioms CRALP.RationalTrace.expTrace
