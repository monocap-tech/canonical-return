import WeilDefect.MasterReturn.LogCarrierNatural
import WeilDefect.MasterReturn.CarrierTranslation

/-! RETURN108-Z: transport the literal physical translation through the
certified continuous carrier equivalence. No signed-form claim. -/
namespace WeilDefect.MasterReturn.LogCarrierTranslation
open MeasureTheory
open scoped FourierTransform
open WeilDefect.MasterRecon
open LogCarrierContinuous LogCarrierNatural
open CarrierTranslation

set_option maxHeartbeats 2000000
set_option maxRecDepth 4096

noncomputable def shiftIV (a b s : ℝ) (hab : a + |s| ≤ b)
    (g : NeutralLogHilbertCarrier a) : NeutralLogHilbertCarrier b :=
  (carrierContinuousEquiv b).symm
    (shiftCarrier a b s hab (carrierContinuousEquiv a g))

theorem shift_natural (a b s : ℝ) (hab : a + |s| ≤ b)
    (g : NeutralLogHilbertCarrier a) :
    carrierContinuousEquiv b (shiftIV a b s hab g) =
      shiftCarrier a b s hab (carrierContinuousEquiv a g) := by
  unfold shiftIV
  exact (carrierContinuousEquiv b).apply_symm_apply _

theorem shiftIV_physical (a b s : ℝ) (hab : a + |s| ≤ b)
    (g : NeutralLogHilbertCarrier a) :
    neutralLogPhysical (shiftIV a b s hab g).val =
      shiftPhysical s (neutralLogPhysical g.val) := by
  rw [← carrierContinuousEquiv_physical b (shiftIV a b s hab g),
    shift_natural, shiftCarrier_physical, carrierContinuousEquiv_physical]

theorem shiftIV_norm (a b s : ℝ) (hab : a + |s| ≤ b)
    (g : NeutralLogHilbertCarrier a) :
    ‖shiftIV a b s hab g‖ = ‖g‖ := by
  have he : ‖shiftIV a b s hab g‖ ^ 2 = ‖g‖ ^ 2 := by
    rw [neutralLogHilbertCarrier_norm_sq (shiftIV a b s hab g),
      neutralLogHilbertCarrier_norm_sq g, shiftIV_physical,
      ← carrierContinuousEquiv_physical a g]
    apply integral_congr_ae
    filter_upwards [shiftPhysical_fourier_norm_ae a b s hab
      (carrierContinuousEquiv a g)] with x hx
    change ‖(𝓕 (shiftPhysical s (V08I.physicalInclusion a
      (carrierContinuousEquiv a g))) : RealComplexL2) x‖ =
      ‖(𝓕 (V08I.physicalInclusion a (carrierContinuousEquiv a g)) :
        RealComplexL2) x‖ at hx
    rw [hx]
  nlinarith [norm_nonneg (shiftIV a b s hab g), norm_nonneg g]

theorem shiftIV_mellin (a b s : ℝ) (hab : a + |s| ≤ b)
    (g : NeutralLogHilbertCarrier a) (z : ℂ) :
    V07A.mellin (neutralLogPhysical (shiftIV a b s hab g).val : ℝ → ℂ) z =
      Complex.exp (z * (s : ℂ)) *
        V07A.mellin (neutralLogPhysical g.val : ℝ → ℂ) z := by
  rw [shiftIV_physical, shiftPhysical_mellin]

theorem shiftIV_target_natural (a b c s : ℝ) (hab : a + |s| ≤ b)
    (hbc : b ≤ c) (g : NeutralLogHilbertCarrier a) :
    ivIncl b c hbc (shiftIV a b s hab g) =
      shiftIV a c s (hab.trans hbc) g := by
  apply LogCarrierBridge.iv_physical_injective c
  change neutralLogPhysical (ivIncl b c hbc (shiftIV a b s hab g)).val =
    neutralLogPhysical (shiftIV a c s (hab.trans hbc) g).val
  rw [ivIncl_physical, shiftIV_physical, shiftIV_physical]

end WeilDefect.MasterReturn.LogCarrierTranslation

#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shiftIV
#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shift_natural
#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shiftIV_physical
#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shiftIV_norm
#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shiftIV_mellin
#print axioms WeilDefect.MasterReturn.LogCarrierTranslation.shiftIV_target_natural
