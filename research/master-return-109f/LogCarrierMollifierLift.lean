import WeilDefect.MasterReturn.LogCarrierCompressionLift
import WeilDefect.MasterReturn.LogCarrierInteriorMollifier

/-! F4c5b: actual smooth compact L2 approximants, logarithmic admission
and exact Fourier multiplier identity. The mollifier sequence is supplied
as an explicit family with radii in the interior collar. -/
namespace WeilDefect.MasterReturn.LogCarrierMollifierLift

open MeasureTheory Filter FourierTransform
open scoped Topology FourierTransform Convolution
open WeilDefect.MasterRecon
open LogCarrierCompressionL2 LogCarrierCompressionLift
open LogCarrierDilationFourier LogCarrierInteriorMollifier
open LogCarrierFormDensityReceivers LogCarrierDilationEnergy
noncomputable section

def smoothWindow (a s : ℝ) (φ : ContDiffBump (0 : ℝ))
    (g : V08I.logCarrier a) : ℝ → ℂ :=
  smoothCompressed φ s (window a g)

theorem smoothWindow_smooth (a s : ℝ) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    ContDiff ℝ ∞ (smoothWindow a s φ g) :=
  smoothCompressed_smooth φ s hs _ (window_integrable a g)

theorem smoothWindow_compact (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    HasCompactSupport (smoothWindow a s φ g) :=
  smoothCompressed_compact φ a s ha hs _ (window_zero a g)

theorem smoothWindow_interior (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) :
    Function.support (smoothWindow a s φ g) ⊆ Set.Ioo (-a) a :=
  smoothCompressed_interior φ a s ha hs hr _ (window_zero a g)

theorem smoothWindow_memLp (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    MemLp (smoothWindow a s φ g) 2 volume :=
  (smoothWindow_smooth a s hs φ g).continuous.memLp_of_hasCompactSupport
    (smoothWindow_compact a s ha hs φ g)

def smoothL2 (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) : V08I.FullL2 :=
  (smoothWindow_memLp a s ha hs φ g).toLp (smoothWindow a s φ g)

theorem smoothL2_ae (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    (smoothL2 a s ha hs φ g : ℝ → ℂ) =ᵐ[volume] smoothWindow a s φ g :=
  (smoothWindow_memLp a s ha hs φ g).coeFn_toLp

theorem smoothL2_support (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) :
    ∀ᵐ x ∂volume, x ∉ Set.Icc (-a) a → smoothL2 a s ha hs φ g x = 0 := by
  filter_upwards [smoothL2_ae a s ha hs φ g] with x hx
  intro hout
  rw [hx]
  by_contra hn
  have hi := smoothWindow_interior a s ha hs φ hr g hn
  exact hout ⟨hi.1.le, hi.2.le⟩

theorem smoothL2_integrable (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    Integrable (smoothL2 a s ha hs φ g : ℝ → ℂ) :=
  ((smoothWindow_smooth a s hs φ g).continuous.integrable_of_hasCompactSupport
    (smoothWindow_compact a s ha hs φ g)).congr
      (smoothL2_ae a s ha hs φ g).symm

theorem smoothL2_fourier_ae (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (φ : ContDiffBump (0 : ℝ)) (g : V08I.logCarrier a) :
    (𝓕 (smoothL2 a s ha hs φ g) : V08I.FullL2) =ᵐ[volume]
      fun u => (𝓕 (bumpComplex φ)) u *
        (𝓕 (compressedL2 a s hs g) : V08I.FullL2) u := by
  filter_upwards [neutralIntegrableL2_fourier_ae (smoothL2 a s ha hs φ g)
    (smoothL2_integrable a s ha hs φ g), compressedL2_fourier_ae a s hs g] with u hu hc
  rw [hu, Real.fourier_congr_ae (smoothL2_ae a s ha hs φ g), hc]
  unfold smoothWindow smoothCompressed
  rw [bump_convolution_fourier φ _ (compress_integrable s hs _ (window_integrable a g)),
    compress_fourier s hs]

theorem smoothL2_log_integrable (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1) (φ : ContDiffBump (0 : ℝ))
    (g : V08I.logCarrier a) :
    Integrable (fun u => V08I.logWeight u *
      ‖(𝓕 (smoothL2 a s ha (by linarith) φ g) : V08I.FullL2) u‖^2) := by
  apply (compressed_log_integrable a s hs hs' g).mono'
    (V08I.logWeight_continuous.aestronglyMeasurable.mul
      (Lp.aestronglyMeasurable (𝓕 (smoothL2 a s ha (by linarith) φ g) :
        V08I.FullL2)).norm.pow 2)
  filter_upwards [smoothL2_fourier_ae a s ha (by linarith) φ g] with u hu
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity), hu, norm_mul]
  have hb := bump_fourier_norm_le φ u
  have hn : ‖(𝓕 (bumpComplex φ)) u‖ *
      ‖(𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2) u‖ ≤
      ‖(𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2) u‖ := by
    simpa only [one_mul] using mul_le_mul_of_nonneg_right hb (norm_nonneg _)
  exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) hn 2)
    (V08I.logWeight_pos u).le

def smoothLift (a s : ℝ) (ha : 0 ≤ a) (hs : (1/2 : ℝ) ≤ s)
    (hs' : s ≤ 1) (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) : V08I.logCarrier a :=
  V08J.liftPhysical a (smoothL2 a s ha (by linarith) φ g)
    (smoothL2_support a s ha (by linarith) φ hr g)
    (smoothL2_log_integrable a s ha hs hs' φ g)

theorem smoothLift_physical (a s : ℝ) (ha : 0 ≤ a) (hs : (1/2 : ℝ) ≤ s)
    (hs' : s ≤ 1) (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) :
    V08I.physicalInclusion a (smoothLift a s ha hs hs' φ hr g) =
      smoothL2 a s ha (by linarith) φ g :=
  V08J.physicalInclusion_lift _ _ _ _

theorem smoothLift_error (a s : ℝ) (ha : 0 ≤ a) (hs : (1/2 : ℝ) ≤ s)
    (hs' : s ≤ 1) (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) :
    ‖smoothLift a s ha hs hs' φ hr g-compression a s ha hs hs' g‖^2 =
      ∫ u, V08I.logWeight u *
      ‖((𝓕 (bumpComplex φ)) u-1) *
        (𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2) u‖^2 := by
  rw [V08I.logCarrier_energy, map_sub, smoothLift_physical, compression_physical,
    FourierTransform.fourier_sub]
  apply integral_congr_ae
  filter_upwards [Lp.coeFn_sub
    (𝓕 (smoothL2 a s ha (by linarith) φ g) : V08I.FullL2)
    (𝓕 (compressedL2 a s (by linarith) g) : V08I.FullL2),
    smoothL2_fourier_ae a s ha (by linarith) φ g] with u hsub hu
  rw [hsub, hu, sub_one_mul]

end
end WeilDefect.MasterReturn.LogCarrierMollifierLift

#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothWindow_smooth
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothWindow_compact
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothWindow_interior
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothWindow_memLp
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothL2_ae
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothL2_support
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothL2_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothL2_fourier_ae
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothL2_log_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothLift_physical
#print axioms WeilDefect.MasterReturn.LogCarrierMollifierLift.smoothLift_error
