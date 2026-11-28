import WeilDefect.MasterReturn.LogCarrierDilationFourier
import Mathlib.Analysis.Calculus.ContDiff.Convolution
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Tactic

/-! F4c3: actual normalized bump, Fourier error and interior support.
No assumption of pointwise smoothness of the original physical vector. -/
namespace WeilDefect.MasterReturn.LogCarrierInteriorMollifier

open MeasureTheory FourierTransform
open scoped FourierTransform Convolution
open LogCarrierDilationFourier
noncomputable section

theorem bumpComplex_smooth (φ : ContDiffBump (0 : ℝ)) :
    ContDiff ℝ ∞ (bumpComplex φ) :=
  Complex.ofRealCLM.contDiff.comp φ.contDiff_normed

theorem bumpComplex_zero (φ : ContDiffBump (0 : ℝ)) (x : ℝ)
    (hx : φ.rOut ≤ |x|) : bumpComplex φ x = 0 := by
  have hz : φ.normed volume x = 0 := by
    by_contra hn
    have hmem : x ∈ Function.support (φ.normed volume) := hn
    rw [φ.support_normed_eq] at hmem
    have hlt : |x| < φ.rOut := by
      simpa [Metric.mem_ball, Real.dist_eq] using hmem
    exact (not_lt_of_ge hx) hlt
  simp [bumpComplex, hz]

theorem bumpComplex_compact (φ : ContDiffBump (0 : ℝ)) :
    HasCompactSupport (bumpComplex φ) := by
  apply HasCompactSupport.intro (K := Set.Icc (-φ.rOut) φ.rOut) isCompact_Icc
  intro x hx
  apply bumpComplex_zero
  exact le_of_not_ge (fun h => hx (abs_le.mp h))

theorem bump_fourier_error (φ : ContDiffBump (0 : ℝ)) (u : ℝ) :
    ‖(𝓕 (bumpComplex φ)) u - 1‖ ≤
      2 * Real.pi * |u| * φ.rOut := by
  let c : ℝ := 2 * Real.pi * |u| * φ.rOut
  let E : ℝ → ℂ := fun x => Complex.exp (Complex.I * ((-2*Real.pi*x*u : ℝ) : ℂ))
  have hE : Continuous E := by dsimp [E]; fun_prop
  have hB : Continuous (bumpComplex φ) := (bumpComplex_smooth φ).continuous
  have hi : Integrable (fun x => E x * bumpComplex φ x) :=
    (hE.mul hB).integrable_of_hasCompactSupport
      ((bumpComplex_compact φ).mul_left E)
  have hmass : (∫ x, bumpComplex φ x) = 1 := by
    simpa only [bumpComplex, Complex.real_smul, mul_one] using
      (φ.integral_normed_smul (μ := volume) (1 : ℂ))
  have he : (𝓕 (bumpComplex φ)) u - 1 =
      ∫ x, E x * bumpComplex φ x - bumpComplex φ x := by
    rw [integral_sub hi (bumpComplex_integrable φ), hmass,
      Real.fourier_real_eq_integral_exp_smul]
    congr 1
    apply integral_congr_ae
    filter_upwards [] with x
    simp only [E, smul_eq_mul, mul_comm Complex.I]
  rw [he]
  apply (norm_integral_le_integral_norm _).trans
  have hbound (x : ℝ) :
      ‖E x * bumpComplex φ x - bumpComplex φ x‖ ≤ c * φ.normed volume x := by
    by_cases hx : φ.rOut ≤ |x|
    · rw [bumpComplex_zero φ x hx]
      simp only [mul_zero, sub_self, norm_zero]
      exact mul_nonneg (by dsimp [c]; positivity) (φ.nonneg_normed x)
    · rw [← sub_one_mul, norm_mul]
      have hphase : ‖E x - 1‖ ≤ |(-2*Real.pi*x*u : ℝ)| := by
        simpa only [E, Real.norm_eq_abs] using
          (Real.norm_exp_I_mul_ofReal_sub_one_le (x := -2*Real.pi*x*u))
      have hphase' : ‖E x - 1‖ ≤ c := by
        apply hphase.trans
        rw [abs_mul, abs_mul, abs_mul, abs_of_nonneg Real.pi_pos.le]
        norm_num
        dsimp [c]
        calc
          2 * Real.pi * |x| * |u| ≤ 2 * Real.pi * φ.rOut * |u| :=
            mul_le_mul_of_nonneg_right
              (mul_le_mul_of_nonneg_left (le_of_lt (lt_of_not_ge hx)) (by positivity))
              (abs_nonneg u)
          _ = _ := by ring
      rw [bumpComplex, Complex.norm_real, Real.norm_eq_abs,
        abs_of_nonneg (φ.nonneg_normed x)]
      exact mul_le_mul_of_nonneg_right hphase' (φ.nonneg_normed x)
  calc
    _ ≤ ∫ x, c * φ.normed volume x :=
      integral_mono (hi.sub (bumpComplex_integrable φ)).norm
        (φ.integrable_normed.const_mul c) hbound
    _ = _ := by rw [integral_const_mul, φ.integral_normed, mul_one]; rfl

def smoothCompressed (φ : ContDiffBump (0 : ℝ)) (s : ℝ) (f : ℝ → ℂ) : ℝ → ℂ :=
  bumpComplex φ ⋆[ContinuousLinearMap.mul ℂ ℂ] compress s f

theorem compress_integrable (s : ℝ) (hs : 0 < s) (f : ℝ → ℂ)
    (hf : Integrable f) : Integrable (compress s f) :=
  hf.comp_mul_left' (inv_ne_zero hs.ne')

theorem compress_zero (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s) (f : ℝ → ℂ)
    (hf : ∀ x, |x| > a → f x = 0) (x : ℝ) (hx : |x| > s*a) :
    compress s f x = 0 := by
  apply hf
  rw [abs_mul, abs_inv, abs_of_pos hs, ← div_eq_inv_mul]
  exact (lt_div_iff₀ hs).mpr (by linarith)

theorem smoothCompressed_smooth (φ : ContDiffBump (0 : ℝ)) (s : ℝ)
    (hs : 0 < s) (f : ℝ → ℂ) (hf : Integrable f) :
    ContDiff ℝ ∞ (smoothCompressed φ s f) := by
  unfold smoothCompressed
  exact (bumpComplex_compact φ).contDiff_convolution_left
    (ContinuousLinearMap.mul ℝ ℂ) (bumpComplex_smooth φ)
    (compress_integrable s hs f hf).locallyIntegrable

theorem smoothCompressed_zero (φ : ContDiffBump (0 : ℝ)) (a s : ℝ)
    (ha : 0 ≤ a) (hs : 0 < s) (f : ℝ → ℂ)
    (hf : ∀ x, |x| > a → f x = 0) (x : ℝ)
    (hx : |x| > s*a + φ.rOut) : smoothCompressed φ s f x = 0 := by
  unfold smoothCompressed
  apply integral_eq_zero_of_ae
  filter_upwards [] with y
  by_cases hy : φ.rOut ≤ |y|
  · simp [bumpComplex_zero φ y hy]
  · have ht : |x| ≤ |x-y| + |y| := by
      calc
        |x| = |(x-y)+y| := by ring_nf
        _ ≤ _ := abs_add_le _ _
    have hz : |x-y| > s*a := by linarith [lt_of_not_ge hy]
    simp [compress_zero a s ha hs f hf _ hz]

end
end WeilDefect.MasterReturn.LogCarrierInteriorMollifier

#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.bumpComplex_smooth
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.bumpComplex_zero
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.bumpComplex_compact
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.bump_fourier_error
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.compress_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.compress_zero
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.smoothCompressed_smooth
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorMollifier.smoothCompressed_zero
