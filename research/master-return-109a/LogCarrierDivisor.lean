import WeilDefect.MasterReturn.LogCarrierTranslation
import WeilDefect.Arithmetic.ActualZetaMultiplicitySymmetry
import WeilDefect.MasterRecon.V10A

/-! RETURN109-A: actual IV raw source evaluations and RECON normalized
orbit samples. Pointwise only: no global lp membership or operator bound. -/
namespace WeilDefect.MasterReturn.LogCarrierDivisor
open MeasureTheory Set InnerProductSpace
open scoped ComplexConjugate InnerProduct ENat
open WeilDefect.MasterRecon
open LogCarrierContinuous

set_option maxHeartbeats 3000000
set_option maxRecDepth 4096

def toPoint (i : V10A.Index) : NeutralActualZetaZeroPoint :=
  ⟨i.val, i.property.1⟩

theorem point_multiplicity (i : V10A.Index) :
    neutralActualZetaMultiplicity (toPoint i) = V10A.multiplicity i := by
  change neutralActualZetaMultiplicity (toPoint i) =
    (analyticOrderAt riemannZeta (toPoint i).val).toNat
  simpa using congrArg ENat.toNat (neutralActualZetaMultiplicity_order (toPoint i))

theorem strip_reflect_closed (ρ : ℂ) (hρ : V09D.stripZero ρ) :
    V09D.stripZero (V10A.reflect ρ) := by
  let q : NeutralActualZetaZeroPoint := ⟨ρ, hρ⟩
  exact (neutralActualZetaZeroPair q).property

theorem reflected_multiplicity (i : V10A.Index) :
    V09D.mult (V10A.reflect (i : ℂ)) = V10A.multiplicity i := by
  change (analyticOrderAt riemannZeta
    (neutralActualZetaZeroPair (toPoint i)).val).toNat =
    (analyticOrderAt riemannZeta (toPoint i).val).toNat
  exact congrArg ENat.toNat ((neutralActualZetaOrder_reflect
    (neutralActualZetaZeroConjugate (toPoint i))).trans
      (neutralActualZetaOrder_conjugate (toPoint i)))
    
theorem window_mellin (a : ℝ) (g : NeutralLogHilbertCarrier a) (z : ℂ) :
    neutralWindowEvaluation a z (neutralLogPhysical g.val) =
      V07A.mellin (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ)
        (Complex.I * z) := by
  rw [carrierContinuousEquiv_physical]
  unfold neutralWindowEvaluation V07A.mellin
  rw [← integral_indicator measurableSet_Icc]
  apply integral_congr_ae
  filter_upwards [neutralLogHilbertCarrier_supported g] with x hx
  by_cases h : x ∈ Icc (-a) a
  · rw [Set.indicator_of_mem h]
    exact mul_comm _ _
  · simp [h, hx h]

private theorem ordinate_exponent (ρ : NeutralActualZetaZeroPoint) :
    Complex.I * neutralActualZetaOrdinate ρ = ρ.val - 1 / 2 := by
  have h := neutralActualZetaOrdinate_source ρ
  linear_combination h

theorem window_right (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) :
    neutralWindowEvaluation a (neutralActualZetaOrdinate (toPoint i))
        (neutralLogPhysical g.val) =
      V07A.mellin (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ)
        ((i : ℂ) - 1 / 2) := by
  rw [window_mellin, ordinate_exponent]
  rfl

theorem window_left (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) :
    neutralWindowEvaluation a (conj (neutralActualZetaOrdinate (toPoint i)))
        (neutralLogPhysical g.val) =
      V07A.mellin (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ)
        (V10A.reflect (i : ℂ) - 1 / 2) := by
  rw [← neutralActualZetaOrdinate_pair, window_mellin, ordinate_exponent]
  rfl

theorem nonfixed_positive (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (hi : V10A.displacement i ≠ 0) :
    (WithLp.ofLp (V10A.normalizedPair
      (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).1 =
      (Real.sqrt (V10A.multiplicity i : ℝ) : ℂ) *
        inner ℂ (neutralLogPositivePairSource a
          (neutralActualZetaOrdinate (toPoint i))) g := by
  rw [V10A.normalizedPair_nonfixed _ _ hi]
  change V07A.pairPlus _ _ = _
  rw [V10A.right_exponent, V10A.left_exponent]
  unfold V07A.weightedMellin V07A.pairPlus
  rw [← window_right a g i, ← window_left a g i, neutralLogPositivePairSource_inner]
  ring

theorem nonfixed_negative (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (hi : V10A.displacement i ≠ 0) :
    (WithLp.ofLp (V10A.normalizedPair
      (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).2 =
      -((Real.sqrt (V10A.multiplicity i : ℝ) : ℂ) *
        inner ℂ (neutralLogNegativePairSource a
          (neutralActualZetaOrdinate (toPoint i))) g) := by
  rw [V10A.normalizedPair_nonfixed _ _ hi]
  change V07A.pairMinus _ _ = _
  rw [V10A.right_exponent, V10A.left_exponent]
  unfold V07A.weightedMellin V07A.pairMinus
  rw [← window_right a g i, ← window_left a g i, neutralLogNegativePairSource_inner]
  ring

theorem fixed_pair (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (hi : V10A.displacement i = 0) :
    WithLp.ofLp (V10A.normalizedPair
      (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i) =
      ((Real.sqrt (V10A.multiplicity i : ℝ) : ℂ) *
        neutralWindowEvaluation a (neutralActualZetaOrdinate (toPoint i))
          (neutralLogPhysical g.val), 0) := by
  rw [V10A.normalizedPair_fixed _ _ hi]
  have hz := V10A.right_exponent i
  rw [hi] at hz
  simp only [Complex.ofReal_zero, zero_add] at hz
  rw [hz]
  unfold V07A.weightedMellin
  rw [← window_right a g i]

end WeilDefect.MasterReturn.LogCarrierDivisor

#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.point_multiplicity
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.strip_reflect_closed
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.reflected_multiplicity
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.window_mellin
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.window_right
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.window_left
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.nonfixed_positive
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.nonfixed_negative
#print axioms WeilDefect.MasterReturn.LogCarrierDivisor.fixed_pair
