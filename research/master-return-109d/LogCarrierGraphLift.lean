import WeilDefect.MasterReturn.LogCarrierGlobalEnergy

/-! RETURN109-D: converse actual graph lift, from literal RECON channel
summability. No universal sampling bound assumed as a theorem. -/
namespace WeilDefect.MasterReturn.LogCarrierGraphLift
open InnerProductSpace
open scoped InnerProduct ComplexConjugate ENNReal
open WeilDefect.MasterRecon
open LogCarrierOrbit LogCarrierGlobalEnergy

noncomputable section

set_option maxHeartbeats 3000000

theorem channel_summability (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hp : Summable (posSq a g)) (hn : Summable (negSq a g)) :
    Summable (fun q : NeutralActualZetaDivisorCoordinate =>
      ‖inner ℂ (neutralLogPositivePairSource a (neutralActualZetaDivisorOrdinate q)) g‖ ^ 2) ∧
    Summable (fun q : NeutralActualZetaDivisorCoordinate =>
      ‖inner ℂ (neutralLogNegativePairSource a (neutralActualZetaDivisorOrdinate q)) g‖ ^ 2) := by
  constructor
  · apply (((positive_summable_iff a g).mpr hp).mul_left (2 : ℝ)).congr
    intro q
    unfold posEnergy
    ring
  · apply (((negative_summable_iff a g).mpr hn).mul_left (2 : ℝ)).congr
    intro q
    unfold negEnergy
    ring

/-- Construct the original graph element, with exactly its raw P/N samples.
Same coordinate construction as the pinned generic square-summable lift;
uses the already compiled graph, avoiding a new import closure. -/
def lift (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hp : Summable (posSq a g)) (hn : Summable (negSq a g)) :
    NeutralActualZetaSourceDomain a := by
  have hs := channel_summability a g hp hn
  let p : NeutralActualZetaGreenCoefficients := ⟨
    fun q => inner ℂ (neutralLogPositivePairSource a (neutralActualZetaDivisorOrdinate q)) g,
    memℓp_gen (by simpa only [ENNReal.toReal_ofNat, Real.rpow_two] using hs.1)⟩
  let n : NeutralActualZetaGreenCoefficients := ⟨
    fun q => inner ℂ (neutralLogNegativePairSource a (neutralActualZetaDivisorOrdinate q)) g,
    memℓp_gen (by simpa only [ENNReal.toReal_ofNat, Real.rpow_two] using hs.2)⟩
  exact ⟨(g,p,n), fun _ => ⟨rfl,rfl⟩⟩

theorem lift_physical (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hp : Summable (posSq a g)) (hn : Summable (negSq a g)) :
    neutralActualZetaSourcePhysical a (lift a g hp hn) = g := rfl

theorem graph_membership_iff (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    (∃ h : NeutralActualZetaSourceDomain a, neutralActualZetaSourcePhysical a h = g) ↔
    Summable (posSq a g) ∧ Summable (negSq a g) := by
  constructor
  · rintro ⟨h,hh⟩
    rw [← hh]
    exact ⟨(graph_positive_energy a h).1,(graph_negative_energy a h).1⟩
  · rintro ⟨hp,hn⟩
    exact ⟨lift a g hp hn,lift_physical a g hp hn⟩

theorem lift_signed (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hp : Summable (posSq a g)) (hn : Summable (negSq a g)) :
    neutralActualZetaSourceQuadratic a (lift a g hp hn) =
      ∑' i, (posSq a g i - negSq a g i) := by
  rw [graph_signed_energy, lift_physical]

theorem lift_coeff_coordinate (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hp : Summable (posSq a g)) (hn : Summable (negSq a g)) (i : V10A.Index) :
    graphCoeff a (lift a g hp hn) i = literalPair a g i := by
  rw [graphCoeff_coordinate, lift_physical]

end
end WeilDefect.MasterReturn.LogCarrierGraphLift

#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.channel_summability
#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.lift
#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.lift_physical
#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.graph_membership_iff
#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.lift_signed
#print axioms WeilDefect.MasterReturn.LogCarrierGraphLift.lift_coeff_coordinate
