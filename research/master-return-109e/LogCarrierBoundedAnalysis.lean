import WeilDefect.MasterReturn.LogCarrierGraphLift
import WeilDefect.MasterRecon.V10B6d
import WeilDefect.External.Zeta23.Theorems.Thm_Zeta23_RvM_zeta_local_zero_count

/-! RETURN109-E: actual local count discharges the existing RECON sampling
receiver, then the admitted D graph lift attaches the original IV graph.
No new estimate, positivity, compactness, or explicit-formula assumption. -/
namespace WeilDefect.MasterReturn.LogCarrierBoundedAnalysis
open InnerProductSpace
open scoped InnerProduct ENNReal
open WeilDefect.MasterRecon
open LogCarrierContinuous LogCarrierGlobalEnergy LogCarrierGraphLift

noncomputable section
set_option maxHeartbeats 3000000
set_option maxRecDepth 4096

theorem count_dictionary (t u : ℝ) :
    V09D.count t u = Zeta23.Ncount t u := rfl

theorem actual_count_bound :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ n : ℤ,
      (V09D.count (n : ℝ) ((n : ℝ)+1) : ℝ) ≤
        A * Real.log (|(n : ℝ)|+3) := by
  obtain ⟨A,hA,h⟩ := Zeta23.RvM.zeta_local_zero_count
  refine ⟨A,by linarith,fun n => ?_⟩
  rw [count_dictionary]
  exact h (n : ℝ)

theorem actual_bounded_analysis (a : ℝ) (ha : 0 ≤ a) :
    ∃ T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index,
      V10B6d.LiteralCoordinates a T := by
  obtain ⟨A,hA,h⟩ := actual_count_bound
  exact V10B6b.bounded_arithmetic_analysis a ha A hA h

theorem literal_channel_summability (a : ℝ)
    (T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index)
    (hT : V10B6d.LiteralCoordinates a T) (g : NeutralLogHilbertCarrier a) :
    Summable (posSq a g) ∧ Summable (negSq a g) := by
  constructor
  · apply (V10B6d.channel_sq_summable true (T (carrierContinuousEquiv a g))).congr
    intro i
    rw [hT (carrierContinuousEquiv a g) i]
    rfl
  · apply (V10B6d.channel_sq_summable false (T (carrierContinuousEquiv a g))).congr
    intro i
    rw [hT (carrierContinuousEquiv a g) i]
    rfl

theorem actual_graph_receiver (a : ℝ) (ha : 0 ≤ a) :
    ∃ T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index,
      V10B6d.LiteralCoordinates a T ∧
      ∀ g : NeutralLogHilbertCarrier a,
        ∃ h : NeutralActualZetaSourceDomain a,
          neutralActualZetaSourcePhysical a h = g ∧
          graphCoeff a h = T (carrierContinuousEquiv a g) ∧
          neutralActualZetaSourceQuadratic a h =
            (inner ℂ (V10B6c.signedGram T (carrierContinuousEquiv a g))
              (carrierContinuousEquiv a g)).re := by
  obtain ⟨T,hT⟩ := actual_bounded_analysis a ha
  refine ⟨T,hT,fun g => ?_⟩
  obtain ⟨hp,hn⟩ := literal_channel_summability a T hT g
  refine ⟨lift a g hp hn,lift_physical a g hp hn,?_,?_⟩
  · apply lp.ext
    funext i
    rw [lift_coeff_coordinate, hT (carrierContinuousEquiv a g) i]
    rfl
  · rw [lift_signed, V10B6c.signedGram_diagonal]
    rw [← V10B6d.signed_channel_tsum (T (carrierContinuousEquiv a g))]
    apply tsum_congr
    intro i
    rw [hT (carrierContinuousEquiv a g) i]
    rfl

end
end WeilDefect.MasterReturn.LogCarrierBoundedAnalysis

#print axioms WeilDefect.MasterReturn.LogCarrierBoundedAnalysis.count_dictionary
#print axioms WeilDefect.MasterReturn.LogCarrierBoundedAnalysis.actual_count_bound
#print axioms WeilDefect.MasterReturn.LogCarrierBoundedAnalysis.actual_bounded_analysis
#print axioms WeilDefect.MasterReturn.LogCarrierBoundedAnalysis.literal_channel_summability
#print axioms WeilDefect.MasterReturn.LogCarrierBoundedAnalysis.actual_graph_receiver
