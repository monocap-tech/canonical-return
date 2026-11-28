import WeilDefect.Arithmetic.ActualZetaCorrelationFourierAgreement
import Mathlib.MeasureTheory.Measure.Haar.NormedSpace
import Mathlib.Analysis.Calculus.BumpFunction.Normed
import Mathlib.Analysis.Fourier.Convolution
import Mathlib.Tactic

/-! RETURN109-F4c2: ordinary Fourier dilation and unit-mass bump
multiplier bounds. Unnormalized physical dilation f(x/s) is enough
for approximation as s tends to one. No density claim. -/
namespace WeilDefect.MasterReturn.LogCarrierDilationFourier

open MeasureTheory FourierTransform
open scoped FourierTransform Convolution
noncomputable section

def compress (s : ℝ) (f : ℝ → ℂ) : ℝ → ℂ :=
  fun x => f (s⁻¹*x)

theorem compress_rawTransform (s : ℝ) (hs : 0 < s)
    (f : ℝ → ℂ) (z : ℂ) :
    neutralRawTransform (compress s f) z =
      (s : ℂ) * neutralRawTransform f ((s : ℂ)*z) := by
  let F : ℝ → ℂ := fun x => f x *
    Complex.exp (Complex.I * ((s : ℂ)*z) * (x : ℂ))
  have h := Measure.integral_comp_inv_mul_left F s
  rw [abs_of_pos hs, Complex.real_smul] at h
  unfold neutralRawTransform compress
  calc
    _ = ∫ x : ℝ, F (s⁻¹*x) := by
      apply integral_congr_ae
      filter_upwards [] with x
      dsimp [F]
      congr 2
      push_cast
      field_simp [hs.ne']
    _ = _ := h

theorem compress_fourier (s : ℝ) (hs : 0 < s)
    (f : ℝ → ℂ) (u : ℝ) :
    (𝓕 (compress s f)) u = (s : ℂ) * (𝓕 f) (s*u) := by
  rw [neutralRawTransform_fourier, compress_rawTransform s hs,
    neutralRawTransform_fourier]
  congr 2
  push_cast
  ring

def bumpComplex (φ : ContDiffBump (0 : ℝ)) : ℝ → ℂ :=
  fun x => (φ.normed volume x : ℂ)

theorem bumpComplex_integrable (φ : ContDiffBump (0 : ℝ)) :
    Integrable (bumpComplex φ) :=
  Complex.ofRealCLM.integrable_comp φ.integrable_normed

theorem bump_fourier_norm_le (φ : ContDiffBump (0 : ℝ)) (u : ℝ) :
    ‖(𝓕 (bumpComplex φ)) u‖ ≤ 1 := by
  rw [Real.fourier_real_eq_integral_exp_smul]
  apply (norm_integral_le_integral_norm _).trans
  have he : (fun x : ℝ => ‖Complex.exp
      ((↑(-2 * Real.pi * x * u) : ℂ) * Complex.I) • bumpComplex φ x‖) =
      φ.normed volume := by
    funext x
    rw [norm_smul, Complex.norm_exp]
    simp [bumpComplex, Complex.mul_re, Complex.norm_real,
      Real.norm_eq_abs, abs_of_nonneg (φ.nonneg_normed x)]
  rw [he, φ.integral_normed]

theorem bump_convolution_fourier (φ : ContDiffBump (0 : ℝ))
    (f : ℝ → ℂ) (hf : Integrable f) (u : ℝ) :
    (𝓕 (bumpComplex φ ⋆[ContinuousLinearMap.mul ℂ ℂ] f)) u =
      (𝓕 (bumpComplex φ)) u * (𝓕 f) u :=
  Real.fourier_mul_convolution_eq (bumpComplex_integrable φ) hf u

end
end WeilDefect.MasterReturn.LogCarrierDilationFourier

#print axioms WeilDefect.MasterReturn.LogCarrierDilationFourier.compress_rawTransform
#print axioms WeilDefect.MasterReturn.LogCarrierDilationFourier.compress_fourier
#print axioms WeilDefect.MasterReturn.LogCarrierDilationFourier.bumpComplex_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierDilationFourier.bump_fourier_norm_le
#print axioms WeilDefect.MasterReturn.LogCarrierDilationFourier.bump_convolution_fourier
