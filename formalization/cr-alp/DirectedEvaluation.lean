import AnalyticGuards
import DirectedInterval
import ElementarySeries
import TrigSeries

namespace CRALP.DirectedEvaluation
open Finset
open scoped BigOperators
open CRALP.DirectedInterval

/-- Transfer a certified polynomial enclosure across its analytic remainder. -/
theorem attach {l u p y e : ℝ} (hp : Contains l u p)
    (he : |y-p| < e) : Contains (l-e) (u+e) y := by
  have h := abs_lt.mp he
  constructor
  · linarith only [hp.1,h.1]
  · linarith only [hp.2,h.2]

theorem logAttach {x l u : ℝ} (hx : 0 < x)
    (hr : |(x-1)/(x+1)| ≤ 7/20)
    (hp : Contains l u (ElementarySeries.logPoly ((x-1)/(x+1)))) :
    Contains (l-ElementarySeries.unit) (u+ElementarySeries.unit) (Real.log x) :=
  attach hp (ElementarySeries.logMöbius x hx hr)

theorem atanAttach {x l u : ℝ} (hx : |x| ≤ 1/3)
    (hp : Contains l u (ElementarySeries.atanPoly 400 x)) :
    Contains (l-ElementarySeries.unit) (u+ElementarySeries.unit) (Real.arctan x) :=
  attach hp (ElementarySeries.atanReduced x hx)

theorem expAttach {x l u : ℝ} (hx : |x| ≤ 3)
    (hp : Contains l u (ElementarySeries.expPoly x)) :
    Contains (l-ElementarySeries.unit) (u+ElementarySeries.unit) (Real.exp x) :=
  attach hp (ElementarySeries.expReduced x hx)

theorem sinAttach {x l u : ℝ} (hx : |x| ≤ 4)
    (hp : Contains l u (TrigSeries.sinPoly 400 x)) :
    Contains (l-TrigSeries.trigUnit) (u+TrigSeries.trigUnit) (Real.sin x) :=
  attach hp (TrigSeries.sinReduced x hx)

theorem cosAttach {x l u : ℝ} (hx : |x| ≤ 4)
    (hp : Contains l u (TrigSeries.cosPoly 400 x)) :
    Contains (l-TrigSeries.trigUnit) (u+TrigSeries.trigUnit) (Real.cos x) :=
  attach hp (TrigSeries.cosReduced x hx)

theorem logHalfAttach {x a b c d : ℝ} (hx : 0 < x)
    (hr : Contains a b (Real.log (x/2))) (h2 : Contains c d (Real.log 2)) :
    Contains (a+c) (b+d) (Real.log x) := by
  rw [ElementarySeries.logHalf x hx]
  exact DirectedInterval.add hr h2

theorem logDoubleAttach {x a b c d : ℝ} (hx : 0 < x)
    (hr : Contains a b (Real.log (2*x))) (h2 : Contains c d (Real.log 2)) :
    Contains (a-d) (b-c) (Real.log x) := by
  rw [ElementarySeries.logDouble x hx]
  exact DirectedInterval.sub hr h2

theorem atanQuarterAttach {x a b c d : ℝ} (hx : 0 ≤ x)
    (hr : Contains a b (Real.arctan ((x-1/4)/(1+x/4))))
    (hq : Contains c d (Real.arctan (1/4))) :
    Contains (a+c) (b+d) (Real.arctan x) := by
  rw [ElementarySeries.atanQuarterReduction x hx]
  exact DirectedInterval.add hr hq

theorem sinPeriodAttach {x a b : ℝ} (n : ℤ)
    (h : Contains a b (Real.sin (x-n*(2*Real.pi)))) :
    Contains a b (Real.sin x) := by
  rwa [ElementarySeries.sinPeriod x n] at h

theorem cosPeriodAttach {x a b : ℝ} (n : ℤ)
    (h : Contains a b (Real.cos (x-n*(2*Real.pi)))) :
    Contains a b (Real.cos x) := by
  rwa [ElementarySeries.cosPeriod x n] at h

theorem machinAttach {a b c d : ℝ}
    (h5 : Contains a b (Real.arctan (1/5)))
    (h239 : Contains c d (Real.arctan (1/239))) :
    Contains (16*a-4*d) (16*b-4*c) Real.pi := by
  rw [ElementarySeries.machin]
  constructor
  · linarith only [h5.1,h239.2]
  · linarith only [h5.2,h239.1]

theorem log2Term (i : ℕ) :
    (1/3 : ℝ)^(2*i+1) / (2*i+1 : ℕ) =
      1 / (((3^(2*i+1)*(2*i+1) : ℕ) : ℝ)) := by
  simp only [Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat, div_pow, one_pow, div_div]

def logDen (i : ℕ) : ℤ := (3^(2*i+1)*(2*i+1) : ℕ)

def logScale : ℤ := 10^340

def logRow (lo hi : ℕ → ℤ) (i : Fin 400) : Prop :=
  lo i.val * logDen i.val ≤ logScale ∧ logScale ≤ hi i.val * logDen i.val

instance (lo hi : ℕ → ℤ) (i : Fin 400) : Decidable (logRow lo hi i) :=
  by unfold logRow; infer_instance

/-- Aggregate exact term bounds and sum endpoints, then apply the proved
    analytic remainder. The actual witness values are supplied separately. -/
theorem log2Trace (lo hi : ℕ → ℤ) (L U : ℤ)
    (hrow : ∀ i : Fin 400, logRow lo hi i)
    (hl : L = 2*(∑ i ∈ range 400, lo i)-1)
    (hu : U = 2*(∑ i ∈ range 400, hi i)+1) :
    Contains ((L : ℝ)/(logScale:ℝ)) ((U : ℝ)/(logScale:ℝ)) (Real.log 2) := by
  have hs : 0 < logScale := by unfold logScale; decide
  have ht (i : ℕ) (hi' : i ∈ range 400) :
      Contains ((lo i : ℝ)/logScale) ((hi i : ℝ)/logScale)
        ((1/3 : ℝ)^(2*i+1)/(2*i+1 : ℕ)) := by
    let j : Fin 400 := ⟨i, mem_range.mp hi'⟩
    have hd : 0 < logDen i := by unfold logDen; positivity
    have h := hrow j
    change lo i * logDen i ≤ logScale ∧ logScale ≤ hi i * logDen i at h
    have hc := DirectedInterval.scaledRatio (lo i) (hi i) 1 (logDen i) logScale hd hs
      (by simpa only [one_mul] using h.1) (by simpa only [one_mul] using h.2)
    rw [log2Term]
    simpa only [logDen, Int.cast_natCast, Int.cast_one] using hc
  have hlo : (∑ i ∈ range 400, (lo i : ℝ)/logScale) ≤
      ∑ i ∈ range 400, (1/3 : ℝ)^(2*i+1)/(2*i+1 : ℕ) :=
    Finset.sum_le_sum (fun i hi' => (ht i hi').1)
  have hhi : (∑ i ∈ range 400, (1/3 : ℝ)^(2*i+1)/(2*i+1 : ℕ)) ≤
      ∑ i ∈ range 400, (hi i : ℝ)/logScale :=
    Finset.sum_le_sum (fun i hi' => (ht i hi').2)
  have hpoly : Contains
      (2*(∑ i ∈ range 400, (lo i : ℝ)/logScale))
      (2*(∑ i ∈ range 400, (hi i : ℝ)/logScale))
      (ElementarySeries.logPoly (1/3)) := by
    unfold ElementarySeries.logPoly
    exact ⟨mul_le_mul_of_nonneg_left hlo (by norm_num),
      mul_le_mul_of_nonneg_left hhi (by norm_num)⟩
  have hr : |(2-1:ℝ)/(2+1)| ≤ 7/20 := by norm_num
  have hp : Contains
      (2*(∑ i ∈ range 400, (lo i : ℝ)/logScale))
      (2*(∑ i ∈ range 400, (hi i : ℝ)/logScale))
      (ElementarySeries.logPoly ((2-1:ℝ)/(2+1))) := by
    rw [show (2-1:ℝ)/(2+1)=1/3 by norm_num]
    exact hpoly
  have hlog := logAttach (x := 2) (by norm_num) hr hp
  have he : ElementarySeries.unit = 1/(logScale:ℝ) := by
    simp only [ElementarySeries.unit, logScale, Nat.cast_pow, Nat.cast_ofNat,
      Int.cast_pow, Int.cast_ofNat]
  have hL : (L : ℝ)/logScale =
      2*(∑ i ∈ range 400, (lo i : ℝ)/logScale)-ElementarySeries.unit := by
    rw [hl,he]
    push_cast
    simp only [div_eq_mul_inv, ← Finset.sum_mul]
    ring
  have hU : (U : ℝ)/logScale =
      2*(∑ i ∈ range 400, (hi i : ℝ)/logScale)+ElementarySeries.unit := by
    rw [hu,he]
    push_cast
    simp only [div_eq_mul_inv, ← Finset.sum_mul]
    ring
  rw [hL,hU]
  exact hlog

end CRALP.DirectedEvaluation
#print axioms CRALP.DirectedEvaluation.attach
#print axioms CRALP.DirectedEvaluation.logAttach
#print axioms CRALP.DirectedEvaluation.atanAttach
#print axioms CRALP.DirectedEvaluation.expAttach
#print axioms CRALP.DirectedEvaluation.sinAttach
#print axioms CRALP.DirectedEvaluation.cosAttach
#print axioms CRALP.DirectedEvaluation.logHalfAttach
#print axioms CRALP.DirectedEvaluation.logDoubleAttach
#print axioms CRALP.DirectedEvaluation.atanQuarterAttach
#print axioms CRALP.DirectedEvaluation.sinPeriodAttach
#print axioms CRALP.DirectedEvaluation.cosPeriodAttach
#print axioms CRALP.DirectedEvaluation.machinAttach
#print axioms CRALP.DirectedEvaluation.log2Term
#print axioms CRALP.DirectedEvaluation.log2Trace
