import AnalyticGuards
import Mathlib

namespace CRALP.ElementarySeries
open Finset Set
open scoped BigOperators
set_option maxRecDepth 20000
set_option maxHeartbeats 0

noncomputable def unit : ℝ := 1 / ((10^340 : ℕ) : ℝ)

noncomputable def logPoly (u : ℝ) : ℝ :=
  2 * ∑ i ∈ range 400, u ^ (2*i+1) / (2*i+1 : ℕ)

noncomputable def atanPoly (n : ℕ) (x : ℝ) : ℝ :=
  ∑ i ∈ range n, (-1 : ℝ)^i * x^(2*i+1) / (2*i+1 : ℕ)

noncomputable def expPoly (x : ℝ) : ℝ :=
  ∑ i ∈ range 400, x^i / (Nat.factorial i : ℝ)

/-- The existing library logarithm bound is slightly looser than the
    historical bound, but the exact original widening still suffices. -/
theorem logLooseGuard : 4 * (7/20 : ℝ)^801 < unit := by
  have h : (4 * 7^801 : ℕ) * 10^340 < 20^801 := by decide
  have h' := AnalyticGuard.ratio_lt_unit _ _ _ (by decide) (by decide) h
  unfold unit
  push_cast at h' ⊢
  rw [div_pow, ← mul_div_assoc]
  exact h'

theorem logReduced (u : ℝ) (hu : |u| ≤ 7/20) :
    |Real.log ((1+u)/(1-u)) - logPoly u| < unit := by
  have h1 : |u| < 1 := lt_of_le_of_lt hu (by norm_num)
  have h := Real.sum_range_sub_log_div_le h1 400
  have hu2 : u^2 ≤ (7/20 : ℝ)^2 := by
    rw [← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) hu 2
  have hd : (1/2 : ℝ) ≤ 1-u^2 := by nlinarith
  have hpow : |u|^801 ≤ (7/20 : ℝ)^801 :=
    pow_le_pow_left₀ (abs_nonneg _) hu 801
  have hb : |u|^801 / (1-u^2) ≤ 2*(7/20 : ℝ)^801 := by
    apply (div_le_iff₀ (by linarith : 0 < 1-u^2)).2
    calc
      |u|^801 ≤ (7/20 : ℝ)^801 := hpow
      _ ≤ 2*(7/20 : ℝ)^801*(1-u^2) := by
        nlinarith only [mul_nonneg (show 0 ≤ (7/20 : ℝ)^801 by positivity)
          (show 0 ≤ 2*(1-u^2)-1 by linarith)]
  have hh : |Real.log ((1+u)/(1-u)) - logPoly u| =
      2 * |1/2 * Real.log ((1+u)/(1-u)) -
        ∑ i ∈ range 400, u^(2*i+1)/(2*i+1 : ℕ)| := by
    rw [← abs_of_nonneg (show (0:ℝ) ≤ 2 by norm_num), ← abs_mul]
    congr 1
    unfold logPoly
    ring
  rw [hh]
  have hbase : |1/2 * Real.log ((1+u)/(1-u)) -
      ∑ i ∈ range 400, u^(2*i+1)/(2*i+1 : ℕ)| ≤
      |u|^801 / (1-u^2) := by
    simpa only [Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat] using h
  calc
    _ ≤ 2 * (2*(7/20 : ℝ)^801) :=
      mul_le_mul_of_nonneg_left (hbase.trans hb) (by norm_num)
    _ = 4*(7/20 : ℝ)^801 := by ring
    _ < unit := logLooseGuard

theorem logMöbius (x : ℝ) (hx : 0 < x)
    (h : |(x-1)/(x+1)| ≤ 7/20) :
    |Real.log x - logPoly ((x-1)/(x+1))| < unit := by
  have he : (1+(x-1)/(x+1))/(1-(x-1)/(x+1)) = x := by
    have hn : x+1 ≠ 0 := by linarith
    have hd : 1-(x-1)/(x+1) ≠ 0 := by
      intro hz
      have hq : (x-1)/(x+1)=1 := by linarith
      have hh := (div_eq_iff hn).1 hq
      linarith
    field_simp [hn,hd]
    <;> ring
  simpa only [he] using logReduced ((x-1)/(x+1)) h

theorem logHalf (x : ℝ) (hx : 0 < x) :
    Real.log x = Real.log (x/2) + Real.log 2 := by
  rw [Real.log_div hx.ne' (by norm_num : (2:ℝ) ≠ 0)]
  ring

theorem logDouble (x : ℝ) (hx : 0 < x) :
    Real.log x = Real.log (2*x) - Real.log 2 := by
  rw [Real.log_mul (by norm_num : (2:ℝ) ≠ 0) hx.ne']
  ring

theorem geomAtan (n : ℕ) (x : ℝ) :
    (1+x^2)*(∑ i ∈ range n, (-1 : ℝ)^i*x^(2*i)) =
      1-(-1 : ℝ)^n*x^(2*n) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, mul_add, ih]
    simp only [pow_succ (-1 : ℝ), Nat.mul_succ, pow_add]
    ring

theorem atanErrorDerivative (n : ℕ) (x : ℝ) :
    HasDerivAt (fun t => Real.arctan t - atanPoly n t)
      ((-1 : ℝ)^n*x^(2*n)/(1+x^2)) x := by
  unfold atanPoly
  have hpoly : HasDerivAt
      (fun t : ℝ => ∑ i ∈ range n, (-1 : ℝ)^i*t^(2*i+1)/(2*i+1 : ℕ))
      (∑ i ∈ range n, (-1 : ℝ)^i*x^(2*i)) x := by
    refine (HasDerivAt.fun_sum fun i hi =>
      ((hasDerivAt_pow (2*i+1) x).const_mul ((-1 : ℝ)^i)).div_const
        ((2*i+1 : ℕ) : ℝ)).congr_deriv ?_
    apply sum_congr rfl
    intro i hi
    simp only [Nat.add_sub_cancel]
    have hne : ((2*i+1 : ℕ) : ℝ) ≠ 0 := by positivity
    field_simp
  refine ((Real.hasDerivAt_arctan x).sub hpoly).congr_deriv ?_
  have hpos : 0 < 1+x^2 := by positivity
  apply (eq_div_iff hpos.ne').2
  linarith only [geomAtan n x, show (1/(1+x^2) : ℝ)*(1+x^2)=1 by field_simp]

theorem atanLooseGuard : (1/3 : ℝ)^801 < unit := by
  have h : (1 : ℕ)*10^340 < 3^801 := by decide
  have h' := AnalyticGuard.ratio_lt_unit _ _ _ (by decide) (by decide) h
  unfold unit
  push_cast at h' ⊢
  simpa only [div_pow, one_pow] using h'

theorem atanReduced (x : ℝ) (hx : |x| ≤ 1/3) :
    |Real.arctan x - atanPoly 400 x| < unit := by
  let F := fun y : ℝ => Real.arctan y - atanPoly 400 y
  let F' := fun y : ℝ => (-1 : ℝ)^400*y^800/(1+y^2)
  have A : ∀ y ∈ Icc (-|x|) |x|, HasDerivAt F (F' y) y := by
    intro y _
    exact atanErrorDerivative 400 y
  have B : ∀ y ∈ Icc (-|x|) |x|, ‖F' y‖ ≤ |x|^800 := by
    intro y hy
    have hpow := pow_le_pow_left₀ (abs_nonneg y) (abs_le.2 hy) 800
    have hd : (1:ℝ) ≤ 1+y^2 := by linarith [sq_nonneg y]
    rw [Real.norm_eq_abs]
    dsimp [F']
    rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_pow, abs_of_pos (show 0 < 1+y^2 by positivity)]
    apply (div_le_iff₀ (by positivity : 0 < 1+y^2)).2
    exact hpow.trans (le_mul_of_one_le_right (by positivity) hd)
  have C := (convex_Icc (-|x|) |x|).norm_image_sub_le_of_norm_hasDerivWithin_le
      (fun y hy => (A y hy).hasDerivWithinAt) B
      (show 0 ∈ Icc (-|x|) |x| by simp)
      (show x ∈ Icc (-|x|) |x| by simp [le_abs_self, neg_le, neg_le_abs x])
  have h0 : F 0 = 0 := by simp [F,atanPoly]
  have hC : |F x| ≤ |x|^801 := by
    rw [pow_succ]
    simpa only [h0, sub_zero, Real.norm_eq_abs] using C
  exact lt_of_le_of_lt (hC.trans (pow_le_pow_left₀ (abs_nonneg _) hx 801)) atanLooseGuard

theorem atanQuarterReduction (x : ℝ) (hx : 0 ≤ x) :
    Real.arctan x = Real.arctan ((x-1/4)/(1+x/4)) + Real.arctan (1/4) := by
  have h := Real.arctan_add (x := x) (y := -(1/4:ℝ)) (by linarith)
  rw [Real.arctan_neg] at h
  have he : (x + -(1/4:ℝ))/(1-x * -(1/4:ℝ)) = (x-1/4)/(1+x/4) := by ring
  rw [he] at h
  linarith

theorem machin :
    Real.pi = 16*Real.arctan (1/5) - 4*Real.arctan (1/239) := by
  have h := Real.four_mul_arctan_inv_5_sub_arctan_inv_239
  simp only [← one_div] at h
  linarith

theorem sinPeriod (x : ℝ) (n : ℤ) :
    Real.sin (x-n*(2*Real.pi)) = Real.sin x :=
  Real.sin_sub_int_mul_two_pi x n

theorem cosPeriod (x : ℝ) (n : ℤ) :
    Real.cos (x-n*(2*Real.pi)) = Real.cos x :=
  Real.cos_sub_int_mul_two_pi x n

theorem expLooseGuard : 3^400 / (Nat.factorial 400 : ℝ) * 2 < unit := by
  have h : (2*3^400 : ℕ)*10^340 < Nat.factorial 400 := by decide
  have h' := AnalyticGuard.ratio_lt_unit _ _ _ (by decide) (by decide) h
  unfold unit
  push_cast at h' ⊢
  calc
    (3:ℝ)^400 / (Nat.factorial 400 : ℝ) * 2 =
        2*3^400/(Nat.factorial 400 : ℝ) := by ring
    _ < _ := h'

theorem expReduced (x : ℝ) (hx : |x| ≤ 3) :
    |Real.exp x-expPoly x| < unit := by
  have hc : ‖(x : ℂ)‖ / ((400 : ℕ).succ : ℝ) ≤ 1/2 := by
    rw [Complex.norm_real]
    change |x| / 401 ≤ 1/2
    linarith
  have h := Complex.exp_bound' (n := 400) hc
  have hr : |Real.exp x - expPoly x| ≤ |x|^400/(Nat.factorial 400 : ℝ)*2 := by
    unfold expPoly
    exact_mod_cast h
  exact lt_of_le_of_lt (hr.trans (by gcongr)) expLooseGuard

end CRALP.ElementarySeries
#print axioms CRALP.ElementarySeries.logLooseGuard
#print axioms CRALP.ElementarySeries.logReduced
#print axioms CRALP.ElementarySeries.logMöbius
#print axioms CRALP.ElementarySeries.logHalf
#print axioms CRALP.ElementarySeries.logDouble
#print axioms CRALP.ElementarySeries.geomAtan
#print axioms CRALP.ElementarySeries.atanErrorDerivative
#print axioms CRALP.ElementarySeries.atanLooseGuard
#print axioms CRALP.ElementarySeries.atanReduced
#print axioms CRALP.ElementarySeries.atanQuarterReduction
#print axioms CRALP.ElementarySeries.machin
#print axioms CRALP.ElementarySeries.sinPeriod
#print axioms CRALP.ElementarySeries.cosPeriod
#print axioms CRALP.ElementarySeries.expLooseGuard
#print axioms CRALP.ElementarySeries.expReduced
