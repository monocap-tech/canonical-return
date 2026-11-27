import WeilDefect.MasterReturn.LogCarrierDivisor

/-! RETURN109-B: exact full divisor-copy indexing and per-orbit signed
energy. No global summability or positivity premise. -/
namespace WeilDefect.MasterReturn.LogCarrierOrbit
open InnerProductSpace
open scoped InnerProduct ComplexConjugate BigOperators
open WeilDefect.MasterRecon
open LogCarrierContinuous LogCarrierDivisor

set_option maxHeartbeats 3000000
set_option maxRecDepth 4096

noncomputable section

abbrev OrbitSlot (i : V10A.Index) :=
  Fin (V10A.multiplicity i) ⊕
    {k : Fin (V10A.multiplicity i) // V10A.displacement i ≠ 0}

abbrev OrbitCopy := Σ i : V10A.Index, OrbitSlot i

def copyRight (i : V10A.Index) (k : Fin (V10A.multiplicity i)) :
    NeutralActualZetaDivisorCoordinate :=
  ⟨toPoint i, ⟨k.val, by rw [point_multiplicity]; exact k.isLt⟩⟩

def copyLeft (i : V10A.Index) (k : Fin (V10A.multiplicity i)) :
    NeutralActualZetaDivisorCoordinate :=
  neutralActualZetaDivisorPair (copyRight i k)

def orbitCopyToDivisor (q : OrbitCopy) : NeutralActualZetaDivisorCoordinate :=
  match q.2 with
  | .inl k => copyRight q.1 k
  | .inr k => copyLeft q.1 k.val

private theorem fixed_of_reflect (i : V10A.Index)
    (h : V10A.reflect (i : ℂ) = (i : ℂ)) : V10A.displacement i = 0 := by
  have he := (V10A.reflect_fixed_iff (i : ℂ)).mp h
  unfold V10A.displacement
  linarith

theorem orbitCopyToDivisor_injective : Function.Injective orbitCopyToDivisor := by
  rintro ⟨i,k⟩ ⟨j,l⟩ h
  rcases k with k | k
  · rcases l with l | l
    · have hp := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.1.val) h
      change (i : ℂ) = (j : ℂ) at hp
      have hij : i = j := Subtype.ext hp
      subst j
      have hk := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.2.val) h
      change k.val = l.val at hk
      have hkl : k = l := Fin.ext hk
      subst l
      rfl
    · have hp := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.1.val) h
      change (i : ℂ) = V10A.reflect (j : ℂ) at hp
      have hij := V10A.orbit_representative_unique i j (Or.inr hp)
      subst j
      exact (l.property (fixed_of_reflect i hp.symm)).elim
  · rcases l with l | l
    · have hp := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.1.val) h
      change V10A.reflect (i : ℂ) = (j : ℂ) at hp
      have hji := V10A.orbit_representative_unique j i (Or.inr hp.symm)
      subst j
      exact (k.property (fixed_of_reflect i hp)).elim
    · have hp := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.1.val) h
      change V10A.reflect (i : ℂ) = V10A.reflect (j : ℂ) at hp
      have he := congrArg V10A.reflect hp
      rw [V10A.reflect_involution, V10A.reflect_involution] at he
      have hij : i = j := Subtype.ext he
      subst j
      have hk := congrArg (fun q : NeutralActualZetaDivisorCoordinate => q.2.val) h
      change k.val.val = l.val.val at hk
      have hkl : k = l := Subtype.ext (Fin.ext hk)
      subst l
      rfl

private theorem copyRight_exists (i : V10A.Index) (q : NeutralActualZetaDivisorCoordinate)
    (h : q.1.val = (i : ℂ)) : ∃ k : Fin (V10A.multiplicity i), copyRight i k = q := by
  have hp : toPoint i = q.1 := Subtype.ext h.symm
  have hm : V10A.multiplicity i = neutralActualZetaMultiplicity q.1 :=
    (point_multiplicity i).symm.trans (congrArg neutralActualZetaMultiplicity hp)
  let k : Fin (V10A.multiplicity i) := ⟨q.2.val, by rw [hm]; exact q.2.isLt⟩
  refine ⟨k, ?_⟩
  apply Sigma.ext hp
  exact (Fin.heq_ext_iff (congrArg neutralActualZetaMultiplicity hp)).mpr rfl

theorem orbitCopyToDivisor_surjective : Function.Surjective orbitCopyToDivisor := by
  intro q
  obtain ⟨i, hi | hi⟩ := V10A.orbit_coverage strip_reflect_closed q.1.val q.1.property
  · obtain ⟨k,hk⟩ := copyRight_exists i q hi
    exact ⟨⟨i,.inl k⟩, hk⟩
  · by_cases hd : V10A.displacement i = 0
    · have hr : (i : ℂ).re = 1 / 2 := by
        unfold V10A.displacement at hd
        linarith
      obtain ⟨k,hk⟩ := copyRight_exists i q
        (hi.trans ((V10A.reflect_fixed_iff (i : ℂ)).mpr hr))
      exact ⟨⟨i,.inl k⟩, hk⟩
    · have hp : (neutralActualZetaDivisorPair q).1.val = (i : ℂ) := by
        change V10A.reflect q.1.val = (i : ℂ)
        rw [hi, V10A.reflect_involution]
      obtain ⟨k,hk⟩ := copyRight_exists i (neutralActualZetaDivisorPair q) hp
      refine ⟨⟨i,.inr ⟨k,hd⟩⟩, ?_⟩
      change neutralActualZetaDivisorPair (copyRight i k) = q
      rw [hk, neutralActualZetaDivisorPair_involutive q]

def orbitCopyEquiv : OrbitCopy ≃ NeutralActualZetaDivisorCoordinate :=
  Equiv.ofBijective orbitCopyToDivisor
    ⟨orbitCopyToDivisor_injective, orbitCopyToDivisor_surjective⟩

/-- The original full-divisor half-factor, before global source-graph assembly. -/
def posEnergy (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (q : NeutralActualZetaDivisorCoordinate) : ℝ :=
  (1 / 2 : ℝ) * ‖inner ℂ
    (neutralLogPositivePairSource a (neutralActualZetaDivisorOrdinate q)) g‖ ^ 2

def negEnergy (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (q : NeutralActualZetaDivisorCoordinate) : ℝ :=
  (1 / 2 : ℝ) * ‖inner ℂ
    (neutralLogNegativePairSource a (neutralActualZetaDivisorOrdinate q)) g‖ ^ 2

private theorem posEnergy_left (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (k : Fin (V10A.multiplicity i)) :
    posEnergy a g (copyLeft i k) = posEnergy a g (copyRight i k) := by
  unfold posEnergy copyLeft
  rw [neutralActualZetaDivisorOrdinate_pair, neutralLogPositivePairSource_conjugate]

private theorem negEnergy_left (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (k : Fin (V10A.multiplicity i)) :
    negEnergy a g (copyLeft i k) = negEnergy a g (copyRight i k) := by
  unfold negEnergy copyLeft
  rw [neutralActualZetaDivisorOrdinate_pair, neutralLogNegativePairSource_conjugate,
    inner_neg_left, norm_neg]

private theorem weighted_norm_sq (m : ℕ) (z : ℂ) :
    ‖(Real.sqrt (m : ℝ) : ℂ) * z‖ ^ 2 = (m : ℝ) * ‖z‖ ^ 2 := by
  rw [norm_mul, mul_pow, Complex.norm_real,
    Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.sq_sqrt (by positivity)]

theorem nonfixed_orbit_energy (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (hi : V10A.displacement i ≠ 0) :
    (∑ k : Fin (V10A.multiplicity i),
      (posEnergy a g (copyRight i k) + posEnergy a g (copyLeft i k))) =
      ‖(WithLp.ofLp (V10A.normalizedPair
        (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).1‖ ^ 2 ∧
    (∑ k : Fin (V10A.multiplicity i),
      (negEnergy a g (copyRight i k) + negEnergy a g (copyLeft i k))) =
      ‖(WithLp.ofLp (V10A.normalizedPair
        (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).2‖ ^ 2 := by
  classical
  constructor
  · rw [nonfixed_positive a g i hi, weighted_norm_sq]
    calc
      _ = ∑ _k : Fin (V10A.multiplicity i),
          ‖inner ℂ (neutralLogPositivePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 := by
        apply Finset.sum_congr rfl
        intro k _
        rw [posEnergy_left]
        unfold posEnergy
        change (1/2 : ℝ) * ‖inner ℂ (neutralLogPositivePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 +
          (1/2 : ℝ) * ‖inner ℂ (neutralLogPositivePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 =
          ‖inner ℂ (neutralLogPositivePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2
        ring
      _ = _ := by simp [nsmul_eq_mul]
  · rw [nonfixed_negative a g i hi, norm_neg, weighted_norm_sq]
    calc
      _ = ∑ _k : Fin (V10A.multiplicity i),
          ‖inner ℂ (neutralLogNegativePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 := by
        apply Finset.sum_congr rfl
        intro k _
        rw [negEnergy_left]
        unfold negEnergy
        change (1/2 : ℝ) * ‖inner ℂ (neutralLogNegativePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 +
          (1/2 : ℝ) * ‖inner ℂ (neutralLogNegativePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2 =
          ‖inner ℂ (neutralLogNegativePairSource a (neutralActualZetaOrdinate (toPoint i))) g‖ ^ 2
        ring
      _ = _ := by simp [nsmul_eq_mul]

private theorem fixed_ordinate (i : V10A.Index) (hi : V10A.displacement i = 0) :
    conj (neutralActualZetaOrdinate (toPoint i)) = neutralActualZetaOrdinate (toPoint i) := by
  apply Complex.ext
  · exact Complex.conj_re _
  · simp only [Complex.conj_im, neutralActualZetaOrdinate_im]
    unfold V10A.displacement at hi
    change -(1 / 2 - (i : ℂ).re) = 1 / 2 - (i : ℂ).re
    linarith

theorem fixed_orbit_energy (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (i : V10A.Index) (hi : V10A.displacement i = 0) :
    (∑ k : Fin (V10A.multiplicity i), posEnergy a g (copyRight i k)) =
      ‖(WithLp.ofLp (V10A.normalizedPair
        (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).1‖ ^ 2 ∧
    (∑ k : Fin (V10A.multiplicity i), negEnergy a g (copyRight i k)) =
      ‖(WithLp.ofLp (V10A.normalizedPair
        (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) i)).2‖ ^ 2 := by
  classical
  let z := neutralActualZetaOrdinate (toPoint i)
  let E := neutralWindowEvaluation a z (neutralLogPhysical g.val)
  have hz : conj z = z := fixed_ordinate i hi
  have hs : (Real.sqrt 2 : ℂ) * (Real.sqrt 2 : ℂ) = 2 := by
    exact_mod_cast Real.mul_self_sqrt (by norm_num : (0:ℝ) ≤ 2)
  have hn : (Real.sqrt 2 : ℂ) ≠ 0 := by
    exact_mod_cast (Real.sqrt_pos.2 (by norm_num : (0:ℝ) < 2)).ne'
  have hscale : (Real.sqrt 2 : ℂ)⁻¹ * 2 = (Real.sqrt 2 : ℂ) := by
    rw [← hs, ← mul_assoc, inv_mul_cancel₀ hn, one_mul]
  have hp : inner ℂ (neutralLogPositivePairSource a z) g = (Real.sqrt 2 : ℂ) * E := by
    rw [neutralLogPositivePairSource_inner, hz]
    change (Real.sqrt 2 : ℂ)⁻¹ * (E + E) = _
    rw [← two_mul, ← mul_assoc, hscale]
  have hm : inner ℂ (neutralLogNegativePairSource a z) g = 0 := by
    rw [neutralLogNegativePairSource_inner, hz, sub_self, mul_zero]
  have hp2 : ‖(Real.sqrt 2 : ℂ) * E‖ ^ 2 = 2 * ‖E‖ ^ 2 := by
    rw [norm_mul, mul_pow, Complex.norm_real,
      Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.sq_sqrt (by norm_num)]
  rw [fixed_pair a g i hi]
  constructor
  · change (∑ k : Fin (V10A.multiplicity i), posEnergy a g (copyRight i k)) =
      ‖(Real.sqrt (V10A.multiplicity i : ℝ) : ℂ) * E‖ ^ 2
    rw [weighted_norm_sq]
    calc
      _ = ∑ _k : Fin (V10A.multiplicity i), ‖E‖ ^ 2 := by
        apply Finset.sum_congr rfl
        intro k _
        change (1/2 : ℝ) * ‖inner ℂ (neutralLogPositivePairSource a z) g‖ ^ 2 = _
        rw [hp, hp2]
        ring
      _ = _ := by simp [nsmul_eq_mul]
  · change (∑ k : Fin (V10A.multiplicity i), negEnergy a g (copyRight i k)) =
      ‖(0 : ℂ)‖ ^ 2
    calc
      _ = ∑ _k : Fin (V10A.multiplicity i), (0 : ℝ) := by
        apply Finset.sum_congr rfl
        intro k _
        change (1/2 : ℝ) * ‖inner ℂ (neutralLogNegativePairSource a z) g‖ ^ 2 = 0
        rw [hm]
        norm_num
      _ = _ := by simp

end
end WeilDefect.MasterReturn.LogCarrierOrbit

#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.copyRight
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.copyLeft
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.orbitCopyToDivisor_injective
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.orbitCopyToDivisor_surjective
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.orbitCopyEquiv
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.nonfixed_orbit_energy
#print axioms WeilDefect.MasterReturn.LogCarrierOrbit.fixed_orbit_energy
