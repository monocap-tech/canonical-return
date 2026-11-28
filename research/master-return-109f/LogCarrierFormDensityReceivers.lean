import WeilDefect.MasterRecon.V08J
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.Analysis.Normed.Lp.SmoothApprox
import Mathlib.Analysis.Calculus.BumpFunction.Normed
import Mathlib.Analysis.Fourier.Convolution
import Mathlib.Tactic

/-! RETURN109-F4c: exact weighted approximation receivers.
No assertion that the physical dilation or mollifier families have already
been constructed or that the full interior core is dense. -/
namespace WeilDefect.MasterReturn.LogCarrierFormDensityReceivers

open MeasureTheory Filter
open scoped FourierTransform Topology
open WeilDefect.MasterRecon
noncomputable section

theorem logWeight_div_le (s u : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) :
    V08I.logWeight (u/s) ≤ V08I.logWeight u + Real.log 2 := by
  have hs0 : 0 < s := by linarith
  have he : |2 * Real.pi * (u/s)| = |2 * Real.pi * u| / s := by
    rw [show 2 * Real.pi * (u/s) = (2 * Real.pi * u)/s by ring,
      abs_div, abs_of_pos hs0]
  have hd : |2 * Real.pi * u| / s ≤ 2 * |2 * Real.pi * u| := by
    apply (div_le_iff₀ hs0).mpr
    nlinarith [mul_nonneg (abs_nonneg (2 * Real.pi * u))
      (sub_nonneg.mpr hs)]
  unfold V08I.logWeight
  rw [he, ← Real.log_mul (by positivity : 2 + |2 * Real.pi * u| ≠ 0)
    (by norm_num : (2 : ℝ) ≠ 0)]
  apply Real.log_le_log (by positivity)
  linarith

theorem logWeight_div_le_twice (s u : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) :
    V08I.logWeight (u/s) ≤ 2 * V08I.logWeight u := by
  linarith [logWeight_div_le s u hs hs', V08I.logWeight_lower u]

theorem logCarrier_sub_energy (a : ℝ) (f g : V08I.logCarrier a) :
    ‖f-g‖^2 = ∫ u, V08I.logWeight u *
      ‖(𝓕 (V08I.physicalInclusion a (f-g)) : V08I.FullL2) u‖^2 :=
  V08I.logCarrier_energy a (f-g)

theorem multiplier_error_bound (z c : ℂ) (hc : ‖c‖ ≤ 1) :
    ‖(c-1)*z‖^2 ≤ 4*‖z‖^2 := by
  have hc' : ‖c-1‖ ≤ 2 := by
    have h := norm_sub_le c (1 : ℂ)
    simp only [norm_one] at h
    linarith
  rw [norm_mul, mul_pow]
  have hsq : ‖c-1‖^2 ≤ 4 := by
    nlinarith [norm_nonneg (c-1),
      mul_nonneg (norm_nonneg (c-1)) (sub_nonneg.mpr hc')]
  exact mul_le_mul_of_nonneg_right hsq (sq_nonneg _)

theorem weighted_multiplier_error_integral_tendsto
    (a : ℝ) (g : V08I.logCarrier a) (m : ℕ → ℝ → ℂ)
    (hm : ∀ n, Measurable (m n))
    (hb : ∀ n u, ‖m n u‖ ≤ 1)
    (ht : ∀ᵐ u ∂volume, Tendsto (fun n => m n u) atTop (𝓝 1)) :
    Tendsto (fun n => ∫ u, V08I.logWeight u *
      ‖(m n u-1) * (𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) u‖^2)
      atTop (𝓝 0) := by
  let F : V08I.FullL2 := 𝓕 (V08I.physicalInclusion a g)
  have hE := V08I.logCarrier_energy_integrable a g
  change Integrable (fun u => V08I.logWeight u * ‖F u‖^2) at hE
  have hD := hE.const_mul 4
  have h := tendsto_integral_of_dominated_convergence
    (μ := volume)
    (F := fun n u => V08I.logWeight u * ‖(m n u-1)*F u‖^2)
    (f := fun _ => (0 : ℝ))
    (fun u => 4*(V08I.logWeight u * ‖F u‖^2))
    (fun n => V08I.logWeight_continuous.aestronglyMeasurable.mul
      (((hm n).aestronglyMeasurable.sub aestronglyMeasurable_const).mul
        (Lp.aestronglyMeasurable F)).norm.pow 2)
    hD
    (fun n => Filter.Eventually.of_forall fun u => by
      rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      have h := mul_le_mul_of_nonneg_left (multiplier_error_bound (F u) (m n u) (hb n u))
        (V08I.logWeight_pos u).le
      calc
        _ ≤ V08I.logWeight u * (4 * ‖F u‖^2) := h
        _ = _ := by ring)
    (by
      filter_upwards [ht] with u hu
      have h0 : Tendsto (fun n => (m n u-1)*F u) atTop (𝓝 0) := by
        simpa using (hu.sub tendsto_const_nhds).mul
          (tendsto_const_nhds (x := F u))
      simpa using (tendsto_const_nhds (x := V08I.logWeight u)).mul
        (h0.norm.pow 2))
  simpa [F] using h

theorem lift_error_energy (a : ℝ) (g : V08I.logCarrier a)
    (f : V08I.FullL2)
    (hs : ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → f x = 0)
    (hf : Integrable (fun u => V08I.logWeight u * ‖(𝓕 f : V08I.FullL2) u‖^2)) :
    ‖V08J.liftPhysical a f hs hf - g‖^2 =
      ∫ u, V08I.logWeight u *
        ‖(𝓕 (f - V08I.physicalInclusion a g) : V08I.FullL2) u‖^2 := by
  rw [V08I.logCarrier_energy, map_sub, V08J.physicalInclusion_lift]

end
end WeilDefect.MasterReturn.LogCarrierFormDensityReceivers

#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.logWeight_div_le
#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.logWeight_div_le_twice
#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.logCarrier_sub_energy
#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.multiplier_error_bound
#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.weighted_multiplier_error_integral_tendsto
#print axioms WeilDefect.MasterReturn.LogCarrierFormDensityReceivers.lift_error_energy
