import WeilDefect.MasterReturn.LogCarrierBridge

/-! RETURN108-X: quantitative continuous equivalence on the actual W carriers.
Generated deterministically by the W/X pinned closure assembler.
No signed-form, sampling, or support-natural assertion is included. -/
namespace WeilDefect.MasterReturn.LogCarrierContinuous
open MeasureTheory
open scoped FourierTransform
open WeilDefect.MasterRecon
open LogWeightComparison LogCarrierBridge

noncomputable def boundRI : ℝ := 1 + |cRI|
noncomputable def boundIR : ℝ := 1 + |cIR| / Real.log 2

theorem boundRI_nonneg : 0 ≤ boundRI := by
  unfold boundRI
  positivity

theorem boundIR_nonneg : 0 ≤ boundIR := by
  unfold boundIR
  have hlog : 0 < Real.log 2 := Real.log_pos (by norm_num)
  positivity

theorem carrierEquiv_norm_sq_le (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    ‖carrierEquiv a g‖ ^ 2 ≤ boundRI * ‖g‖ ^ 2 := by
  have hi := neutralLogPhysical_energy g.val
  have hr := iv_energy_to_recon _ hi
  rw [V08I.logCarrier_energy, carrierEquiv_physical, neutralLogHilbertCarrier_norm_sq]
  calc
    _ ≤ ∫ x, boundRI * (logarithmicFourierWeight x *
      ‖(𝓕 (neutralLogPhysical g.val) : RealComplexL2) x‖ ^ 2) := by
        apply integral_mono hr (hi.const_mul boundRI)
        intro x
        change V08I.logWeight x * (‖(𝓕 (neutralLogPhysical g.val) : RealComplexL2) x‖ ^ 2) ≤
          boundRI * (logarithmicFourierWeight x * (‖(𝓕 (neutralLogPhysical g.val) : RealComplexL2) x‖ ^ 2))
        rw [← mul_assoc]
        exact mul_le_mul_of_nonneg_right (recon_le_scaled_iv x) (sq_nonneg _)
    _ = _ := by rw [integral_const_mul]

theorem carrierEquiv_symm_norm_sq_le (a : ℝ) (g : V08I.logCarrier a) :
    ‖(carrierEquiv a).symm g‖ ^ 2 ≤ boundIR * ‖g‖ ^ 2 := by
  have hp : neutralLogPhysical ((carrierEquiv a).symm g).val =
      V08I.physicalInclusion a g := toIV_physical a g
  have hr := V08I.logCarrier_energy_integrable a g
  have hi := recon_energy_to_iv _ hr
  rw [neutralLogHilbertCarrier_norm_sq, hp, V08I.logCarrier_energy]
  calc
    _ ≤ ∫ x, boundIR * (V08I.logWeight x *
      ‖(𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) x‖ ^ 2) := by
        apply integral_mono hi (hr.const_mul boundIR)
        intro x
        change logarithmicFourierWeight x * (‖(𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) x‖ ^ 2) ≤
          boundIR * (V08I.logWeight x * (‖(𝓕 (V08I.physicalInclusion a g) : V08I.FullL2) x‖ ^ 2))
        rw [← mul_assoc]
        exact mul_le_mul_of_nonneg_right (iv_le_scaled_recon x) (sq_nonneg _)
    _ = _ := by rw [integral_const_mul]

private theorem nonneg_le_of_pow_two_le {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y)
    (h : x ^ 2 ≤ y ^ 2) : x ≤ y := by
  by_contra! hn
  have hpos : 0 < (x-y)*(x+y) :=
    mul_pos (sub_pos.mpr hn) (add_pos_of_pos_of_nonneg (lt_of_le_of_lt hy hn) hy)
  nlinarith only [hpos, h]

theorem carrierEquiv_norm_le (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    ‖carrierEquiv a g‖ ≤ Real.sqrt boundRI * ‖g‖ := by
  apply nonneg_le_of_pow_two_le (norm_nonneg _) (mul_nonneg (Real.sqrt_nonneg _) (norm_nonneg _))
  rw [mul_pow, Real.sq_sqrt boundRI_nonneg]
  exact carrierEquiv_norm_sq_le a g

theorem carrierEquiv_symm_norm_le (a : ℝ) (g : V08I.logCarrier a) :
    ‖(carrierEquiv a).symm g‖ ≤ Real.sqrt boundIR * ‖g‖ := by
  apply nonneg_le_of_pow_two_le (norm_nonneg _) (mul_nonneg (Real.sqrt_nonneg _) (norm_nonneg _))
  rw [mul_pow, Real.sq_sqrt boundIR_nonneg]
  exact carrierEquiv_symm_norm_sq_le a g

noncomputable def carrierContinuousEquiv (a : ℝ) :
    NeutralLogHilbertCarrier a ≃L[ℂ] V08I.logCarrier a where
  toLinearEquiv := carrierEquiv a
  continuous_toFun :=
    ((carrierEquiv a).toLinearMap.mkContinuous (Real.sqrt boundRI)
      (carrierEquiv_norm_le a)).continuous
  continuous_invFun :=
    ((carrierEquiv a).symm.toLinearMap.mkContinuous (Real.sqrt boundIR)
      (carrierEquiv_symm_norm_le a)).continuous

theorem carrierContinuousEquiv_physical (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    V08I.physicalInclusion a (carrierContinuousEquiv a g) = neutralLogPhysical g.val :=
  carrierEquiv_physical a g

end WeilDefect.MasterReturn.LogCarrierContinuous

#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.boundRI_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.boundIR_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierEquiv_norm_sq_le
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierEquiv_symm_norm_sq_le
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierEquiv_norm_le
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierEquiv_symm_norm_le
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierContinuousEquiv
#print axioms WeilDefect.MasterReturn.LogCarrierContinuous.carrierContinuousEquiv_physical
