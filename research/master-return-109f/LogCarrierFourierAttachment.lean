import WeilDefect.Arithmetic.ActualZetaArchimedeanTransport
import WeilDefect.MasterRecon.V08I
import Mathlib.Tactic

/-! RETURN109-F4a: remove the Green restriction from the literal
correlation Fourier identity. No graph, smoothness, or density premise. -/
namespace WeilDefect.MasterReturn.LogCarrierFourierAttachment

open MeasureTheory FourierTransform
open scoped FourierTransform ComplexConjugate
noncomputable section

theorem supported_window_ae (a : ℝ) (f : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0) :
    neutralWindowRepresentative a f =ᵐ[volume] f := by
  filter_upwards [hf] with x hx
  by_cases hm : x ∈ Set.Icc (-a) a
  · simp [neutralWindowRepresentative, hm]
  · simp [neutralWindowRepresentative, hm, hx hm]

theorem supported_integrable (a : ℝ) (f : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0) :
    Integrable (f : ℝ → ℂ) := by
  have hi : Integrable (neutralWindowRepresentative a f) := by
    simpa using neutralWindowRepresentative_twist_integrable a f 0
  exact hi.congr (supported_window_ae a f hf)

theorem supported_fourier_window (a : ℝ) (f : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0) (ξ : ℝ) :
    (𝓕 (f : ℝ → ℂ)) ξ =
      neutralWindowEvaluation a ((-2 * Real.pi * ξ : ℝ) : ℂ) f := by
  calc
    _ = (𝓕 (neutralWindowRepresentative a f)) ξ :=
      Real.fourier_congr_ae (supported_window_ae a f hf).symm ξ
    _ = _ := neutralWindowRepresentative_fourier a f ξ

theorem correlation_fourier_ae (a : ℝ) (f g : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0)
    (hg : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → g x = 0) :
    (𝓕 (neutralWindowCorrelation a f g)) =ᵐ[volume]
      fun ξ => conj ((𝓕 f : RealComplexL2) ξ) * (𝓕 g : RealComplexL2) ξ := by
  filter_upwards [neutralIntegrableL2_fourier_ae f (supported_integrable a f hf),
    neutralIntegrableL2_fourier_ae g (supported_integrable a g hg)] with ξ hξf hξg
  rw [neutralRawTransform_fourier, neutralWindowCorrelation_rawTransform]
  simp only [Complex.conj_ofReal]
  rw [hξf, hξg, supported_fourier_window a f hf, supported_fourier_window a g hg]

theorem literal_arch_integrand_ae (a : ℝ) (f g : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0)
    (hg : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → g x = 0) :
    (fun ξ : ℝ => Zeta23.paperFT (neutralWindowCorrelation a f g)
      ((-2 * Real.pi * ξ : ℝ) : ℂ) *
      (Zeta23.EF.gammaBracket (-2 * Real.pi * ξ) : ℂ)) =ᵐ[volume]
    fun ξ => (Zeta23.EF.gammaBracket (2 * Real.pi * ξ) : ℂ) *
      (conj ((𝓕 f : RealComplexL2) ξ) * (𝓕 g : RealComplexL2) ξ) := by
  filter_upwards [correlation_fourier_ae a f g hf hg] with ξ hξ
  rw [← neutralRawTransform_eq_paperFT, ← neutralRawTransform_fourier, hξ]
  rw [show -2 * Real.pi * ξ = -(2 * Real.pi * ξ) by ring,
    neutralActualZetaGammaBracket_even]
  exact mul_comm _ _

theorem literal_arch_spectral (a : ℝ) (f g : RealComplexL2)
    (hf : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0)
    (hg : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → g x = 0) :
    (1 / (2 * Real.pi) : ℂ) *
      (∫ r : ℝ, Zeta23.paperFT (neutralWindowCorrelation a f g) (r : ℂ) *
        (Zeta23.EF.gammaBracket r : ℂ)) =
    ∫ ξ : ℝ, (Zeta23.EF.gammaBracket (2 * Real.pi * ξ) : ℂ) *
      (conj ((𝓕 f : RealComplexL2) ξ) * (𝓕 g : RealComplexL2) ξ) := by
  let F : ℝ → ℂ := fun r => Zeta23.paperFT (neutralWindowCorrelation a f g) (r : ℂ) *
      (Zeta23.EF.gammaBracket r : ℂ)
  have h := Measure.integral_comp_mul_left F (-2 * Real.pi)
  have habs : |(-2 * Real.pi)⁻¹| = 1 / (2 * Real.pi) := by
    rw [show -2 * Real.pi = -(2 * Real.pi) by ring,
      inv_neg, abs_neg, abs_of_pos (inv_pos.mpr (by positivity))]
    simp only [one_div]
  rw [habs, Complex.real_smul] at h
  have he : (∫ ξ : ℝ, F ((-2 * Real.pi) * ξ)) = _ :=
    integral_congr_ae (literal_arch_integrand_ae a f g hf hg)
  have hc : (1 / (2 * Real.pi) : ℂ) = ((1 / (2 * Real.pi) : ℝ) : ℂ) := by
    push_cast
    rfl
  rw [hc]
  exact h.symm.trans he

theorem gammaBracket_measurable : Measurable Zeta23.EF.gammaBracket :=
  neutralActualZetaGammaBracket_continuous.measurable

end
end WeilDefect.MasterReturn.LogCarrierFourierAttachment

#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.supported_window_ae
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.supported_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.supported_fourier_window
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.correlation_fourier_ae
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.literal_arch_integrand_ae
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.literal_arch_spectral
#print axioms WeilDefect.MasterReturn.LogCarrierFourierAttachment.gammaBracket_measurable
