import WeilDefect.MasterReturn.LogCarrierDilationConvergence
import WeilDefect.MasterReturn.LogCarrierFormDensityReceivers

/-! F4c5a: bind the strong weighted dilation receiver to the actual
supported spectral carrier. Source candidate, not kernel admission. -/
namespace WeilDefect.MasterReturn.LogCarrierCompressionLift

open MeasureTheory Filter FourierTransform
open scoped Topology FourierTransform
open WeilDefect.MasterRecon
open LogCarrierCompressionL2 LogCarrierDilationEnergy
open LogCarrierDilationConvergence LogCarrierFormDensityReceivers
noncomputable section

theorem compression_support_old (a s : ℝ) (ha : 0 ≤ a)
    (hs : 0 < s) (hs' : s ≤ 1) (g : V08I.logCarrier a) :
    ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → compressedL2 a s hs g x = 0 := by
  filter_upwards [compressedL2_support a s ha hs g] with x hx
  intro hout
  apply hx
  intro hin
  apply hout
  have h := mul_le_mul_of_nonneg_right hs' ha
  exact ⟨by linarith [hin.1], by linarith [hin.2]⟩

def compression (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) (g : V08I.logCarrier a) :
    V08I.logCarrier a :=
  V08J.liftPhysical a (compressedL2 a s (by linarith) g)
    (compression_support_old a s ha (by linarith) hs' g)
    (compressed_log_integrable a s hs hs' g)

theorem compression_physical (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) (g : V08I.logCarrier a) :
    V08I.physicalInclusion a (compression a s ha hs hs' g) =
      compressedL2 a s (by linarith) g :=
  V08J.physicalInclusion_lift _ _ _ _

theorem compression_error (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) (g : V08I.logCarrier a) :
    ‖compression a s ha hs hs' g-g‖^2 =
      ∫ u, err s (𝓕 (window a g)) u := by
  rw [V08I.logCarrier_energy, map_sub, compression_physical, FourierTransform.fourier_sub]
  apply integral_congr_ae
  filter_upwards [Lp.coeFn_sub
    (𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2)
    (𝓕 (V08I.physicalInclusion a g) : V08I.FullL2),
    compressedL2_fourier_ae a s (by linarith) g, window_fourier_ae a g] with u hsub hc hg
  rw [hsub, hc, hg]
  rfl

theorem compression_tendsto (a : ℝ) (ha : 0 ≤ a)
    (g : V08I.logCarrier a) (s : ℕ → ℝ)
    (hs : ∀ n, (1/2 : ℝ) ≤ s n) (hs' : ∀ n, s n ≤ 1)
    (ht : Tendsto s atTop (𝓝 1)) :
    Tendsto (fun n => compression a (s n) ha (hs n) (hs' n) g)
      atTop (𝓝 g) := by
  let B := ∫ x, ‖window a g x‖
  have hB : 0 ≤ B := integral_nonneg fun x => norm_nonneg _
  have hb (u : ℝ) : ‖(𝓕 (window a g)) u‖ ≤ B :=
    VectorFourier.norm_fourierIntegral_le_integral_norm
      Real.fourierChar volume (innerₗ ℝ) (window a g) u
  have he := strong_weighted_dilation (𝓕 (window a g))
    (window_fourier_continuous a g) B hB hb (window_energy_integrable a g)
    s hs hs' ht
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩ := Metric.tendsto_atTop.mp he (ε^2) (by positivity)
  refine ⟨N, fun n hn => ?_⟩
  have h := hN n hn
  rw [← compression_error a (s n) ha (hs n) (hs' n) g,
    Real.dist_eq, sub_zero, abs_of_nonneg (sq_nonneg _)] at h
  rw [dist_eq_norm]
  by_contra hbad
  have hh : ε ≤ ‖compression a (s n) ha (hs n) (hs' n) g-g‖ :=
    le_of_not_gt hbad
  nlinarith [mul_nonneg (sub_nonneg.mpr hh)
    (show 0 ≤ ‖compression a (s n) ha (hs n) (hs' n) g-g‖+ε by positivity)]

end
end WeilDefect.MasterReturn.LogCarrierCompressionLift

#print axioms WeilDefect.MasterReturn.LogCarrierCompressionLift.compression_support_old
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionLift.compression_physical
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionLift.compression_error
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionLift.compression_tendsto
