import WeilDefect.MasterReturn.LogCarrierCompressionL2
import WeilDefect.MasterReturn.LogCarrierFormDensityReceivers
import Mathlib.MeasureTheory.Integral.Bochner.Set

/-! F4c4b: weighted dilation energy and uniform tail receiver.
This is a source candidate, not an admission of strong convergence. -/
namespace WeilDefect.MasterReturn.LogCarrierDilationEnergy

open MeasureTheory FourierTransform Filter
open scoped FourierTransform Topology
open WeilDefect.MasterRecon
open LogCarrierCompressionL2 LogCarrierFormDensityReceivers
noncomputable section

def energy (F : ℝ → ℂ) (u : ℝ) : ℝ :=
  V08I.logWeight u * ‖F u‖^2

def dilatedEnergy (s : ℝ) (F : ℝ → ℂ) (u : ℝ) : ℝ :=
  V08I.logWeight u * ‖(s : ℂ) * F (s*u)‖^2

theorem energy_nonneg (F : ℝ → ℂ) (u : ℝ) : 0 ≤ energy F u :=
  mul_nonneg (V08I.logWeight_pos u).le (sq_nonneg _)

theorem dilatedEnergy_nonneg (s : ℝ) (F : ℝ → ℂ) (u : ℝ) :
    0 ≤ dilatedEnergy s F u :=
  mul_nonneg (V08I.logWeight_pos u).le (sq_nonneg _)

theorem dilatedEnergy_pointwise (s : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (F : ℝ → ℂ) (u : ℝ) :
    dilatedEnergy s F u ≤ 2*s^2 * energy F (s*u) := by
  have hs0 : 0 < s := by linarith
  have hw := logWeight_div_le_twice s (s*u) hs hs'
  rw [show s*u/s = u by field_simp [hs0.ne']] at hw
  unfold energy dilatedEnergy
  rw [norm_mul, mul_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hs0]
  calc
    _ ≤ (2 * V08I.logWeight (s*u)) * (s^2 * ‖F (s*u)‖^2) :=
      mul_le_mul_of_nonneg_right hw (by positivity)
    _ = _ := by ring

theorem dilatedEnergy_integrable (s : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (F : ℝ → ℂ) (hF : Continuous F) (hE : Integrable (energy F)) :
    Integrable (dilatedEnergy s F) := by
  have hs0 : 0 < s := by linarith
  have h := (hE.comp_mul_left' hs0.ne').const_mul (2*s^2)
  apply h.mono' (show Continuous (dilatedEnergy s F) from by
    unfold dilatedEnergy
    fun_prop).aestronglyMeasurable
  filter_upwards [] with u
  rw [Real.norm_eq_abs, abs_of_nonneg (dilatedEnergy_nonneg s F u)]
  exact dilatedEnergy_pointwise s hs hs' F u

theorem dilatedEnergy_budget (s : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (F : ℝ → ℂ) (hE : Integrable (energy F)) :
    (∫ u, dilatedEnergy s F u) ≤ 2*s * ∫ u, energy F u := by
  have hs0 : 0 < s := by linarith
  calc
    _ ≤ ∫ u, 2*s^2 * energy F (s*u) :=
      integral_mono_of_nonneg (Filter.Eventually.of_forall (dilatedEnergy_nonneg s F))
        ((hE.comp_mul_left' hs0.ne').const_mul (2*s^2))
        (Filter.Eventually.of_forall (dilatedEnergy_pointwise s hs hs' F))
    _ = _ := by
      rw [integral_const_mul, Measure.integral_comp_mul_left,
        abs_of_pos (inv_pos.mpr hs0), smul_eq_mul]
      field_simp [hs0.ne']

def tail (R : ℝ) : Set ℝ := {u | R < |u|}

theorem tail_measurable (R : ℝ) : MeasurableSet (tail R) :=
  measurableSet_lt measurable_const continuous_abs.measurable

theorem dilatedEnergy_tail (s R : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (hR : 0 ≤ R) (F : ℝ → ℂ) (hE : Integrable (energy F)) :
    (∫ u in tail R, dilatedEnergy s F u) ≤
      2*s * ∫ u in tail (R/2), energy F u := by
  have hs0 : 0 < s := by linarith
  let H := (tail (R/2)).indicator (energy F)
  have hH : Integrable H := hE.integrableOn.integrable_indicator (tail_measurable _)
  rw [← integral_indicator (tail_measurable R),
    ← integral_indicator (tail_measurable (R/2))]
  calc
    _ ≤ ∫ u, 2*s^2 * H (s*u) := by
      apply integral_mono_of_nonneg
      · filter_upwards [] with u
        by_cases hu : u ∈ tail R <;> simp [Set.indicator, hu, dilatedEnergy_nonneg]
      · exact (hH.comp_mul_left' hs0.ne').const_mul _
      · filter_upwards [] with u
        by_cases hu : u ∈ tail R
        · have hu' : R < |u| := hu
          have ht : s*u ∈ tail (R/2) := by
            change R/2 < |s*u|
            rw [abs_mul, abs_of_pos hs0]
            nlinarith [mul_nonneg (sub_nonneg.mpr hs) (abs_nonneg u)]
          simp only [Set.indicator_of_mem hu, H, Set.indicator_of_mem ht]
          exact dilatedEnergy_pointwise s hs hs' F u
        · rw [Set.indicator_of_notMem hu]
          dsimp [H]
          by_cases ht : s*u ∈ tail (R/2)
          · rw [Set.indicator_of_mem ht]
            exact mul_nonneg (by positivity) (energy_nonneg F _)
          · rw [Set.indicator_of_notMem ht, mul_zero]
    _ = _ := by
      change _ = 2*s * ∫ u, H u
      rw [integral_const_mul, Measure.integral_comp_mul_left,
        abs_of_pos (inv_pos.mpr hs0), smul_eq_mul]
      field_simp [hs0.ne']

theorem window_fourier_continuous (a : ℝ) (g : V08I.logCarrier a) :
    Continuous (𝓕 (window a g)) :=
  VectorFourier.fourierIntegral_continuous (L := innerₗ ℝ)
    Real.continuous_fourierChar continuous_inner (window_integrable a g)

theorem window_energy_integrable (a : ℝ) (g : V08I.logCarrier a) :
    Integrable (energy (𝓕 (window a g))) := by
  apply (V08I.logCarrier_energy_integrable a g).congr
  filter_upwards [window_fourier_ae a g] with u hu
  simp only [energy, hu]

theorem compressed_log_integrable (a s : ℝ) (hs : (1/2 : ℝ) ≤ s)
    (hs' : s ≤ 1) (g : V08I.logCarrier a) :
    Integrable (fun u => V08I.logWeight u *
      ‖(𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2) u‖^2) := by
  apply (dilatedEnergy_integrable s hs hs' _ (window_fourier_continuous a g)
    (window_energy_integrable a g)).congr
  filter_upwards [compressedL2_fourier_ae a s (by linarith) g] with u hu
  simp only [dilatedEnergy, hu]

end
end WeilDefect.MasterReturn.LogCarrierDilationEnergy

#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.energy_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.dilatedEnergy_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.dilatedEnergy_pointwise
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.dilatedEnergy_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.dilatedEnergy_budget
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.tail_measurable
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.dilatedEnergy_tail
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.window_fourier_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.window_energy_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierDilationEnergy.compressed_log_integrable
