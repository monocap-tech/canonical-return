import WeilDefect.MasterReturn.LogCarrierInteriorMollifier
import WeilDefect.MasterReturn.LogCarrierFourierAttachment
import WeilDefect.MasterRecon.V08J
import Mathlib.MeasureTheory.Function.L2Space

/-! F4c4a: actual physical L2 compression of the literal window
representative. No strong convergence assumption or density claim. -/
namespace WeilDefect.MasterReturn.LogCarrierCompressionL2

open MeasureTheory FourierTransform
open scoped FourierTransform
open WeilDefect.MasterRecon
open LogCarrierDilationFourier LogCarrierInteriorMollifier
open LogCarrierFourierAttachment
noncomputable section

theorem compress_memLp (s : ℝ) (hs : 0 < s) (f : ℝ → ℂ)
    (hm : Measurable f) (hf : MemLp f 2 volume) :
    MemLp (compress s f) 2 volume := by
  have hc : Measurable (compress s f) := hm.comp (by fun_prop)
  apply (memLp_two_iff_integrable_sq_norm hc.aestronglyMeasurable).mpr
  exact ((memLp_two_iff_integrable_sq_norm hf.aestronglyMeasurable).mp hf).comp_mul_left'
    (inv_ne_zero hs.ne')

def window (a : ℝ) (g : V08I.logCarrier a) : ℝ → ℂ :=
  neutralWindowRepresentative a (V08I.physicalInclusion a g)

theorem window_measurable (a : ℝ) (g : V08I.logCarrier a) :
    Measurable (window a g) :=
  (V08I.physicalInclusion a g).stronglyMeasurable.measurable.indicator measurableSet_Icc

theorem window_memLp (a : ℝ) (g : V08I.logCarrier a) :
    MemLp (window a g) 2 volume :=
  (Lp.memLp (V08I.physicalInclusion a g)).ae_eq
    (supported_window_ae a _ (V08I.physicalInclusion_support a g)).symm

theorem window_integrable (a : ℝ) (g : V08I.logCarrier a) :
    Integrable (window a g) := by
  simpa [window] using neutralWindowRepresentative_twist_integrable a
    (V08I.physicalInclusion a g) 0

theorem window_zero (a : ℝ) (g : V08I.logCarrier a) (x : ℝ)
    (hx : a < |x|) : window a g x = 0 := by
  have h : x ∉ Set.Icc (-a) a := by
    intro hmem
    exact (not_le_of_gt hx) (abs_le.mpr hmem)
  simp [window, neutralWindowRepresentative, h]

def compressedL2 (a s : ℝ) (hs : 0 < s) (g : V08I.logCarrier a) :
    V08I.FullL2 :=
  (compress_memLp s hs _ (window_measurable a g) (window_memLp a g)).toLp
    (compress s (window a g))

theorem compressedL2_ae (a s : ℝ) (hs : 0 < s) (g : V08I.logCarrier a) :
    (compressedL2 a s hs g : ℝ → ℂ) =ᵐ[volume] compress s (window a g) :=
  (compress_memLp s hs _ (window_measurable a g) (window_memLp a g)).coeFn_toLp

theorem compressedL2_integrable (a s : ℝ) (hs : 0 < s) (g : V08I.logCarrier a) :
    Integrable (compressedL2 a s hs g : ℝ → ℂ) :=
  (compress_integrable s hs _ (window_integrable a g)).congr
    (compressedL2_ae a s hs g).symm

theorem compressedL2_support (a s : ℝ) (ha : 0 ≤ a) (hs : 0 < s)
    (g : V08I.logCarrier a) :
    ∀ᵐ x ∂volume, x ∉ Set.Icc (-(s*a)) (s*a) →
      compressedL2 a s hs g x = 0 := by
  filter_upwards [compressedL2_ae a s hs g] with x hx
  intro hout
  rw [hx]
  exact compress_zero a s ha hs _ (window_zero a g) x
    (lt_of_not_ge (fun h => hout (abs_le.mp h)))

theorem compressedL2_fourier_ae (a s : ℝ) (hs : 0 < s) (g : V08I.logCarrier a) :
    (𝓕 (compressedL2 a s hs g) : V08I.FullL2) =ᵐ[volume]
      fun u => (s : ℂ) * (𝓕 (window a g)) (s*u) := by
  filter_upwards [neutralIntegrableL2_fourier_ae (compressedL2 a s hs g)
    (compressedL2_integrable a s hs g)] with u hu
  rw [hu, Real.fourier_congr_ae (compressedL2_ae a s hs g),
    compress_fourier s hs]

theorem window_fourier_ae (a : ℝ) (g : V08I.logCarrier a) :
    (𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) =ᵐ[volume]
      (𝓕 (window a g)) := by
  filter_upwards [neutralIntegrableL2_fourier_ae (V08I.physicalInclusion a g)
    (supported_integrable a _ (V08I.physicalInclusion_support a g))] with u hu
  rw [hu]
  exact Real.fourier_congr_ae
    (supported_window_ae a _ (V08I.physicalInclusion_support a g)).symm u

end
end WeilDefect.MasterReturn.LogCarrierCompressionL2

#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.compress_memLp
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.window_measurable
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.window_memLp
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.window_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.window_zero
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.compressedL2_ae
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.compressedL2_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.compressedL2_support
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.compressedL2_fourier_ae
#print axioms WeilDefect.MasterReturn.LogCarrierCompressionL2.window_fourier_ae
