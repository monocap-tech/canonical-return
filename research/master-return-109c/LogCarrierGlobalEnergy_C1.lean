import WeilDefect.MasterReturn.LogCarrierOrbit
import WeilDefect.Arithmetic.ActualZetaSelectedBackground
import Mathlib.Topology.Algebra.InfiniteSum.Real

/-! RETURN109-C: global summation on the authentic complete source graph.
Each nonnegative channel is summed separately; no cancellation of divergent
series, no bounded physical analysis or positivity claim. -/
namespace WeilDefect.MasterReturn.LogCarrierGlobalEnergy
open InnerProductSpace
open scoped InnerProduct ComplexConjugate ENNReal BigOperators
open WeilDefect.MasterRecon
open LogCarrierContinuous LogCarrierDivisor LogCarrierOrbit

set_option maxHeartbeats 3000000
set_option maxRecDepth 4096

noncomputable section

def literalPair (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) : V07C.Pair :=
  V10A.normalizedPair
    (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i

def posSq (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) : ℝ :=
  ‖(WithLp.ofLp (literalPair a g i)).1‖ ^ 2

def negSq (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) : ℝ :=
  ‖(WithLp.ofLp (literalPair a g i)).2‖ ^ 2

private theorem fiber_sum (e : NeutralActualZetaDivisorCoordinate → ℝ)
    (i : V10A.Index) :
    (∑' k : OrbitSlot i, e (orbitCopyEquiv ⟨i,k⟩)) =
      if V10A.displacement i = 0 then
        ∑ k : Fin (V10A.multiplicity i), e (copyRight i k)
      else ∑ k : Fin (V10A.multiplicity i),
        (e (copyRight i k) + e (copyLeft i k)) := by
  classical
  rw [tsum_fintype, Fintype.sum_sum_type]
  change ((∑ k : Fin (V10A.multiplicity i), e (copyRight i k)) +
    ∑ k : {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0},
      e (copyLeft i k.val)) = _
  by_cases hd : V10A.displacement i = 0
  · letI : IsEmpty {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0} :=
      ⟨fun k => k.property hd⟩
    simp [hd]
  · let E : Fin (V10A.multiplicity i) ≃
        {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0} :=
      { toFun := fun k => ⟨k,hd⟩
        invFun := Subtype.val
        left_inv := fun _ => rfl
        right_inv := fun k => Subtype.ext rfl }
    have hs : (∑ k : {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0},
        e (copyLeft i k.val)) =
        ∑ k : Fin (V10A.multiplicity i), e (copyLeft i k) := by
      have ht := (E.tsum_eq (fun k => e (copyLeft i k.val))).symm
      change (∑' k : {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0},
        e (copyLeft i k.val)) =
        ∑' k : Fin (V10A.multiplicity i), e (copyLeft i k) at ht
      simpa only [tsum_fintype] using ht
    rw [if_neg hd, hs, Finset.sum_add_distrib]

theorem fiber_positive (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) :
    (∑' k : OrbitSlot i, posEnergy a g (orbitCopyEquiv ⟨i,k⟩)) = posSq a g i := by
  rw [fiber_sum]
  by_cases hd : V10A.displacement i = 0
  · rw [if_pos hd]
    exact (fixed_orbit_energy a g i hd).1
  · rw [if_neg hd]
    exact (nonfixed_orbit_energy a g i hd).1

theorem fiber_negative (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) :
    (∑' k : OrbitSlot i, negEnergy a g (orbitCopyEquiv ⟨i,k⟩)) = negSq a g i := by
  rw [fiber_sum]
  by_cases hd : V10A.displacement i = 0
  · rw [if_pos hd]
    exact (fixed_orbit_energy a g i hd).2
  · rw [if_neg hd]
    exact (nonfixed_orbit_energy a g i hd).2

private theorem orbit_summable_iff
    (e : NeutralActualZetaDivisorCoordinate → ℝ) (v : V10A.Index → ℝ)
    (he : ∀ q, 0 ≤ e q)
    (hv : ∀ i, (∑' k : OrbitSlot i, e (orbitCopyEquiv ⟨i,k⟩)) = v i) :
    Summable e ↔ Summable v := by
  classical
  let F : OrbitCopy → ℝ := e ∘ orbitCopyEquiv
  have hF : Summable F ↔ Summable e := orbitCopyEquiv.summable_iff
  have hσ := summable_sigma_of_nonneg (f := F) (fun q => he (orbitCopyEquiv q))
  constructor
  · intro h
    have hh := (hσ.mp (hF.mpr h)).2
    exact hh.congr hv
  · intro h
    apply hF.mp
    apply hσ.mpr
    refine ⟨fun i => ?_, ?_⟩
    · exact Summable.of_finite
    · exact h.congr (fun i => (hv i).symm)

private theorem orbit_tsum
    (e : NeutralActualZetaDivisorCoordinate → ℝ) (v : V10A.Index → ℝ)
    (he : Summable e)
    (hv : ∀ i, (∑' k : OrbitSlot i, e (orbitCopyEquiv ⟨i,k⟩)) = v i) :
    (∑' q, e q) = ∑' i, v i := by
  let F : OrbitCopy → ℝ := e ∘ orbitCopyEquiv
  have hF : Summable F := orbitCopyEquiv.summable_iff.mpr he
  calc
    (∑' q, e q) = ∑' q : OrbitCopy, F q := (orbitCopyEquiv.tsum_eq e).symm
    _ = ∑' i, ∑' k : OrbitSlot i, F ⟨i,k⟩ := hF.tsum_sigma
    _ = ∑' i, v i := tsum_congr hv

theorem positive_summable_iff (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    Summable (posEnergy a g) ↔ Summable (posSq a g) :=
  orbit_summable_iff _ _ (fun q => by unfold posEnergy; positivity) (fiber_positive a g)

theorem negative_summable_iff (a : ℝ) (g : NeutralLogHilbertCarrier a) :
    Summable (negEnergy a g) ↔ Summable (negSq a g) :=
  orbit_summable_iff _ _ (fun q => by unfold negEnergy; positivity) (fiber_negative a g)

private theorem graph_pos_copy (a : ℝ) (h : NeutralActualZetaSourceDomain a)
    (q : NeutralActualZetaDivisorCoordinate) :
    ‖neutralActualZetaNormalizedPositiveAnalysis a h q‖ ^ 2 =
      posEnergy a (neutralActualZetaSourcePhysical a h) q := by
  simp only [neutralActualZetaNormalizedPositiveAnalysis, ContinuousLinearMap.smul_apply,
    lp.coeFn_smul, Pi.smul_apply, norm_smul, mul_pow,
    neutralActualZetaSourceNormalization_norm_sq]
  rw [(neutralActualZetaSourceAnalysis_sampling a h q).1] <;> rfl

private theorem graph_neg_copy (a : ℝ) (h : NeutralActualZetaSourceDomain a)
    (q : NeutralActualZetaDivisorCoordinate) :
    ‖neutralActualZetaNormalizedNegativeAnalysis a h q‖ ^ 2 =
      negEnergy a (neutralActualZetaSourcePhysical a h) q := by
  simp only [neutralActualZetaNormalizedNegativeAnalysis, ContinuousLinearMap.smul_apply,
    lp.coeFn_smul, Pi.smul_apply, norm_smul, mul_pow,
    neutralActualZetaSourceNormalization_norm_sq]
  rw [(neutralActualZetaSourceAnalysis_sampling a h q).2] <;> rfl

theorem graph_positive_energy (a : ℝ) (h : NeutralActualZetaSourceDomain a) :
    Summable (posSq a (neutralActualZetaSourcePhysical a h)) ∧
    ‖neutralActualZetaNormalizedPositiveAnalysis a h‖ ^ 2 =
      ∑' i, posSq a (neutralActualZetaSourcePhysical a h) i := by
  have hs := lp.hasSum_norm (by norm_num : 0 < (2 : ℝ≥0∞).toReal)
    (neutralActualZetaNormalizedPositiveAnalysis a h)
  simp only [ENNReal.toReal_ofNat, Real.rpow_two] at hs
  have he : Summable (posEnergy a (neutralActualZetaSourcePhysical a h)) :=
    hs.summable.congr (graph_pos_copy a h)
  refine ⟨(positive_summable_iff _ _).mp he, ?_⟩
  calc
    _ = ∑' q, ‖neutralActualZetaNormalizedPositiveAnalysis a h q‖ ^ 2 := hs.tsum_eq.symm
    _ = ∑' q, posEnergy a (neutralActualZetaSourcePhysical a h) q :=
      tsum_congr (graph_pos_copy a h)
    _ = _ := orbit_tsum _ _ he (fiber_positive _ _)

theorem graph_negative_energy (a : ℝ) (h : NeutralActualZetaSourceDomain a) :
    Summable (negSq a (neutralActualZetaSourcePhysical a h)) ∧
    ‖neutralActualZetaNormalizedNegativeAnalysis a h‖ ^ 2 =
      ∑' i, negSq a (neutralActualZetaSourcePhysical a h) i := by
  have hs := lp.hasSum_norm (by norm_num : 0 < (2 : ℝ≥0∞).toReal)
    (neutralActualZetaNormalizedNegativeAnalysis a h)
  simp only [ENNReal.toReal_ofNat, Real.rpow_two] at hs
  have he : Summable (negEnergy a (neutralActualZetaSourcePhysical a h)) :=
    hs.summable.congr (graph_neg_copy a h)
  refine ⟨(negative_summable_iff _ _).mp he, ?_⟩
  calc
    _ = ∑' q, ‖neutralActualZetaNormalizedNegativeAnalysis a h q‖ ^ 2 := hs.tsum_eq.symm
    _ = ∑' q, negEnergy a (neutralActualZetaSourcePhysical a h) q :=
      tsum_congr (graph_neg_copy a h)
    _ = _ := orbit_tsum _ _ he (fiber_negative _ _)

theorem graph_signed_energy (a : ℝ) (h : NeutralActualZetaSourceDomain a) :
    neutralActualZetaSourceQuadratic a h =
      ∑' i, (posSq a (neutralActualZetaSourcePhysical a h) i -
        negSq a (neutralActualZetaSourcePhysical a h) i) := by
  have hp := graph_positive_energy a h
  have hn := graph_negative_energy a h
  unfold neutralActualZetaSourceQuadratic
  rw [hp.2, hn.2, ← hp.1.tsum_sub hn.1]

private theorem pair_norm_sq (a : ℝ) (g : NeutralLogHilbertCarrier a) (i : V10A.Index) :
    ‖literalPair a g i‖ ^ 2 = posSq a g i + negSq a g i := by
  rw [WithLp.prod_norm_eq_of_L2]
  change (Real.sqrt (posSq a g i + negSq a g i)) ^ 2 = _
  exact Real.sq_sqrt (by unfold posSq negSq; positivity)

def graphCoeff (a : ℝ) (h : NeutralActualZetaSourceDomain a) :
    V07C.Coeff V10A.Index := by
  let g := neutralActualZetaSourcePhysical a h
  have hs : Summable (fun i => ‖literalPair a g i‖ ^ 2) :=
    ((graph_positive_energy a h).1.add (graph_negative_energy a h).1).congr
      (fun i => (pair_norm_sq a g i).symm)
  exact ⟨literalPair a g, memℓp_gen (by
    simpa only [ENNReal.toReal_ofNat, Real.rpow_two] using hs)⟩

theorem graphCoeff_coordinate (a : ℝ) (h : NeutralActualZetaSourceDomain a) (i : V10A.Index) :
    graphCoeff a h i = literalPair a (neutralActualZetaSourcePhysical a h) i := rfl

theorem graphCoeff_norm_sq (a : ℝ) (h : NeutralActualZetaSourceDomain a) :
    ‖graphCoeff a h‖ ^ 2 =
      ‖neutralActualZetaNormalizedPositiveAnalysis a h‖ ^ 2 +
      ‖neutralActualZetaNormalizedNegativeAnalysis a h‖ ^ 2 := by
  have hs := (lp.hasSum_norm (by norm_num : 0 < (2 : ℝ≥0∞).toReal)
    (graphCoeff a h)).tsum_eq
  simp only [ENNReal.toReal_ofNat, Real.rpow_two] at hs
  have hp := graph_positive_energy a h
  have hn := graph_negative_energy a h
  calc
    _ = ∑' i, ‖graphCoeff a h i‖ ^ 2 := hs.symm
    _ = ∑' i, (posSq a (neutralActualZetaSourcePhysical a h) i +
      negSq a (neutralActualZetaSourcePhysical a h) i) := by
        apply tsum_congr
        intro i
        rw [graphCoeff_coordinate, pair_norm_sq]
    _ = _ := by rw [hp.1.tsum_add hn.1, ← hp.2, ← hn.2]

end
end WeilDefect.MasterReturn.LogCarrierGlobalEnergy

#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.fiber_positive
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.fiber_negative
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.positive_summable_iff
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.negative_summable_iff
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graph_positive_energy
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graph_negative_energy
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graph_signed_energy
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graphCoeff
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graphCoeff_coordinate
#print axioms WeilDefect.MasterReturn.LogCarrierGlobalEnergy.graphCoeff_norm_sq
