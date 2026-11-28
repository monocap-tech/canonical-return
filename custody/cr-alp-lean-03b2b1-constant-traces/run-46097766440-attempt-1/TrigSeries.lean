import AnalyticGuards
import Mathlib

namespace CRALP.TrigSeries
open Finset Set
open scoped BigOperators
set_option maxRecDepth 20000
set_option maxHeartbeats 0

noncomputable def trigUnit : ℝ := 1 / ((10^340 : ℕ) : ℝ)

noncomputable def derivPoly (f : ℝ → ℝ) (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range n, iteratedDeriv k f 0 * x^k / (Nat.factorial k : ℝ)

noncomputable def sinPoly (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range n, (-1 : ℝ)^k*x^(2*k+1)/(Nat.factorial (2*k+1) : ℝ)

noncomputable def cosPoly (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range n, (-1 : ℝ)^k*x^(2*k)/(Nat.factorial (2*k) : ℝ)

theorem sinPoly_eq (n : ℕ) (x : ℝ) : derivPoly Real.sin (2*n) x = sinPoly n x := by
  induction n with
  | zero => simp [derivPoly,sinPoly]
  | succ n ih =>
    unfold derivPoly sinPoly at *
    rw [show 2*(n+1)=2*n+1+1 by omega, sum_range_succ, sum_range_succ, ih, sum_range_succ]
    simp [Real.iteratedDeriv_even_sin,Real.iteratedDeriv_odd_sin]

theorem cosPoly_eq (n : ℕ) (x : ℝ) : derivPoly Real.cos (2*n) x = cosPoly n x := by
  induction n with
  | zero => simp [derivPoly,cosPoly]
  | succ n ih =>
    unfold derivPoly cosPoly at *
    rw [show 2*(n+1)=2*n+1+1 by omega, sum_range_succ, sum_range_succ, ih, sum_range_succ]
    simp [Real.iteratedDeriv_even_cos,Real.iteratedDeriv_odd_cos]

theorem taylor800 (f : ℝ → ℝ) (hf : ContDiff ℝ 800 f)
    (hb : ∀ y : ℝ, |iteratedDeriv 800 f y| ≤ 1) (x : ℝ) (hx : x ≠ 0) :
    |f x - derivPoly f 800 x| ≤ |x|^800 / (Nat.factorial 800 : ℝ) := by
  have hne : (0:ℝ) ≠ x := Ne.symm hx
  have hf' : ContDiffOn ℝ ((799 : ℕ)+1) f (uIcc 0 x) := hf.contDiffOn
  obtain ⟨y,hy,hr⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 799) hne hf'
  have hpoly : taylorWithinEval f 799 (uIcc 0 x) 0 x = derivPoly f 800 x := by
    rw [taylor_within_apply]
    unfold derivPoly
    apply sum_congr rfl
    intro k hk
    have he : iteratedDerivWithin k f (uIcc 0 x) 0 = iteratedDeriv k f 0 :=
      iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_uIcc hne)
        (hf.of_le (by exact_mod_cast (show k ≤ 800 by
          have hh := mem_range.mp hk
          omega))).contDiffAt (by simp)
    rw [he]
    simp only [sub_zero,smul_eq_mul]
    ring
  rw [hpoly] at hr
  rw [hr,sub_zero,abs_div,abs_mul,abs_pow,abs_of_nonneg (by positivity : 0 ≤ (Nat.factorial 800 : ℝ))]
  exact div_le_div_of_nonneg_right
    (by simpa only [one_mul] using mul_le_mul_of_nonneg_right (hb y) (by positivity : 0 ≤ |x|^800))
    (by positivity)

theorem trigGuard : (4:ℝ)^800 / (Nat.factorial 800 : ℝ) < trigUnit := by
  have h := AnalyticGuard.sincos_400_remainder_real
  unfold trigUnit
  push_cast at h ⊢
  exact h

theorem sinReduced (x : ℝ) (hx : |x| ≤ 4) :
    |Real.sin x-sinPoly 400 x| < trigUnit := by
  by_cases hz : x=0
  · subst x
    simp [sinPoly,trigUnit]
  · have h := taylor800 Real.sin Real.contDiff_sin (Real.abs_iteratedDeriv_sin_le_one 800) x hz
    rw [show (800:ℕ)=2*400 by decide,sinPoly_eq] at h
    exact lt_of_le_of_lt (h.trans (by gcongr)) trigGuard

theorem cosReduced (x : ℝ) (hx : |x| ≤ 4) :
    |Real.cos x-cosPoly 400 x| < trigUnit := by
  by_cases hz : x=0
  · subst x
    have h : cosPoly 400 0 = 1 := by
      unfold cosPoly
      rw [Finset.sum_eq_single 0]
      · simp
      · intro b hb hb0
        have hn : 2*b ≠ 0 := by omega
        simp [zero_pow hn]
      · simp
    rw [h,Real.cos_zero,sub_self,abs_zero]
    unfold trigUnit
    positivity
  · have h := taylor800 Real.cos Real.contDiff_cos (Real.abs_iteratedDeriv_cos_le_one 800) x hz
    rw [show (800:ℕ)=2*400 by decide,cosPoly_eq] at h
    exact lt_of_le_of_lt (h.trans (by gcongr)) trigGuard

end CRALP.TrigSeries
#print axioms CRALP.TrigSeries.sinPoly_eq
#print axioms CRALP.TrigSeries.cosPoly_eq
#print axioms CRALP.TrigSeries.taylor800
#print axioms CRALP.TrigSeries.trigGuard
#print axioms CRALP.TrigSeries.sinReduced
#print axioms CRALP.TrigSeries.cosReduced
