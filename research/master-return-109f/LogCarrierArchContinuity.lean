import WeilDefect.MasterReturn.LogCarrierFourierAttachment
import WeilDefect.MasterReturn.LogCarrierGammaGrowth
import Mathlib.Analysis.InnerProductSpace.Continuous

/-! RETURN109-F4b: actual archimedean diagonal on the exact weighted
spectral carrier, as a bounded multiplier pairing. No zero-side definition. -/
namespace WeilDefect.MasterReturn.LogCarrierArchContinuity

open MeasureTheory InnerProductSpace
open scoped InnerProduct FourierTransform ComplexConjugate
open WeilDefect.MasterRecon
open LogCarrierFourierAttachment LogCarrierGammaGrowth
noncomputable section

def gammaRatio (u : ℝ) : ℂ :=
  ((Zeta23.EF.gammaBracket (2 * Real.pi * u) / V08I.logWeight u : ℝ) : ℂ)

theorem gammaRatio_measurable : Measurable gammaRatio := by
  unfold gammaRatio
  exact Complex.continuous_ofReal.measurable.comp
    ((gammaBracket_measurable.comp (by fun_prop)).div
      V08I.logWeight_continuous.measurable)

theorem gammaRatio_bound (K : ℝ)
    (h : ∀ u : ℝ, |Zeta23.EF.gammaBracket (2 * Real.pi * u)| ≤
      K * V08I.logWeight u) (u : ℝ) :
    ‖gammaRatio u‖ ≤ K := by
  rw [gammaRatio, Complex.norm_real, Real.norm_eq_abs, abs_div,
    abs_of_pos (V08I.logWeight_pos u)]
  exact (div_le_iff₀ (V08I.logWeight_pos u)).mpr (h u)

def gammaOperator (K : ℝ)
    (h : ∀ u : ℝ, |Zeta23.EF.gammaBracket (2 * Real.pi * u)| ≤
      K * V08I.logWeight u) : V08I.FullL2 →L[ℂ] V08I.FullL2 :=
  V08I.multiplier gammaRatio gammaRatio_measurable K (gammaRatio_bound K h)

def literalArch (a : ℝ) (g : V08I.logCarrier a) : ℂ :=
  (1 / (2 * Real.pi) : ℂ) *
    ∫ r : ℝ, Zeta23.paperFT
      (neutralWindowCorrelation a (V08I.physicalInclusion a g)
        (V08I.physicalInclusion a g)) (r : ℂ) *
        (Zeta23.EF.gammaBracket r : ℂ)

theorem literalArch_integrable (a : ℝ) (g : V08I.logCarrier a)
    (K : ℝ) (hK : 0 ≤ K)
    (h : ∀ u : ℝ, |Zeta23.EF.gammaBracket (2 * Real.pi * u)| ≤
      K * V08I.logWeight u) :
    Integrable (fun r : ℝ => Zeta23.paperFT
      (neutralWindowCorrelation a (V08I.physicalInclusion a g)
        (V08I.physicalInclusion a g)) (r : ℂ) *
        (Zeta23.EF.gammaBracket r : ℂ)) := by
  let f := V08I.physicalInclusion a g
  let H : ℝ → ℂ := fun u => (Zeta23.EF.gammaBracket (2 * Real.pi * u) : ℂ) *
    (conj ((𝓕 f : V08I.FullL2) u) * (𝓕 f : V08I.FullL2) u)
  have hm : AEStronglyMeasurable H volume := by
    apply ((Complex.continuous_ofReal.measurable.comp
      (gammaBracket_measurable.comp (by fun_prop))).aestronglyMeasurable).mul
    exact (Lp.aestronglyMeasurable (𝓕 f : V08I.FullL2)).star.mul
      (Lp.aestronglyMeasurable (𝓕 f : V08I.FullL2))
  have hi : Integrable H := by
    apply ((V08I.logCarrier_energy_integrable a g).const_mul K).mono' hm
    filter_upwards [] with u
    dsimp [H, f]
    rw [norm_mul, norm_mul, Complex.norm_conj, Complex.norm_real, Real.norm_eq_abs]
    simpa only [mul_assoc, pow_two] using
      mul_le_mul_of_nonneg_right (h u)
        (mul_nonneg (norm_nonneg _) (norm_nonneg _))
  have hi' := hi.congr (literal_arch_integrand_ae a f f
    (V08I.physicalInclusion_support a g) (V08I.physicalInclusion_support a g)).symm
  exact (integrable_comp_mul_left_iff _ (show -2 * Real.pi ≠ 0 by
    exact mul_ne_zero (by norm_num) Real.pi_ne_zero)).mp hi'

theorem literalArch_pairing (a : ℝ) (g : V08I.logCarrier a)
    (K : ℝ) (h : ∀ u : ℝ, |Zeta23.EF.gammaBracket (2 * Real.pi * u)| ≤
      K * V08I.logWeight u) :
    literalArch a g = inner ℂ (g : V08I.FullL2)
      (gammaOperator K h (g : V08I.FullL2)) := by
  unfold literalArch
  rw [literal_arch_spectral a _ _ (V08I.physicalInclusion_support a g)
    (V08I.physicalInclusion_support a g), L2.inner_def]
  apply integral_congr_ae
  filter_upwards [V08I.logCarrier_energy_ae a g,
    V08I.multiplier_ae gammaRatio gammaRatio_measurable K
      (gammaRatio_bound K h) (g : V08I.FullL2)] with u he hm
  change V08I.logWeight u * ‖(𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) u‖^2 =
    ‖(g : V08I.FullL2) u‖^2 at he
  change gammaOperator K h (g : V08I.FullL2) u =
    gammaRatio u * (g : V08I.FullL2) u at hm
  let F : ℂ := (𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) u
  let z : ℂ := (g : V08I.FullL2) u
  let G : ℝ := Zeta23.EF.gammaBracket (2 * Real.pi * u)
  let w : ℝ := V08I.logWeight u
  have hw : w ≠ 0 := (V08I.logWeight_pos u).ne'
  have henergy : w * ‖F‖^2 = ‖z‖^2 := he
  have hreal : G * ‖F‖^2 = (G / w) * ‖z‖^2 := by
    rw [← henergy]
    field_simp [hw] <;> ring
  have hF : conj F * F = ((‖F‖^2 : ℝ) : ℂ) := by
    rw [Complex.conj_mul'] <;> push_cast
  have hz : conj z * z = ((‖z‖^2 : ℝ) : ℂ) := by
    rw [Complex.conj_mul']
    push_cast
  rw [hm, RCLike.inner_apply']
  change (G : ℂ) * (conj F * F) = conj z * ((↑(G / w) : ℂ) * z)
  calc
    _ = (G : ℂ) * ((‖F‖^2 : ℝ) : ℂ) := by rw [hF]
    _ = ((G / w : ℝ) : ℂ) * ((‖z‖^2 : ℝ) : ℂ) := by
      exact_mod_cast hreal
    _ = conj z * (((G / w : ℝ) : ℂ) * z) := by rw [← hz]; ring

theorem literalArch_continuous (a : ℝ) : Continuous (literalArch a) := by
  obtain ⟨K,hK,h⟩ := gamma_logWeight_bound
  have he : literalArch a = fun g : V08I.logCarrier a =>
      inner ℂ (g : V08I.FullL2) (gammaOperator K h (g : V08I.FullL2)) := by
    funext g
    exact literalArch_pairing a g K h
  rw [he]
  exact continuous_subtype_val.inner
    ((gammaOperator K h).continuous.comp continuous_subtype_val)

end
end WeilDefect.MasterReturn.LogCarrierArchContinuity

#print axioms WeilDefect.MasterReturn.LogCarrierArchContinuity.gammaRatio_measurable
#print axioms WeilDefect.MasterReturn.LogCarrierArchContinuity.gammaRatio_bound
#print axioms WeilDefect.MasterReturn.LogCarrierArchContinuity.literalArch_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierArchContinuity.literalArch_pairing
#print axioms WeilDefect.MasterReturn.LogCarrierArchContinuity.literalArch_continuous
