import WeilDefect.Arithmetic.ActualZetaCorrelationPoles
import WeilDefect.Arithmetic.ActualZetaIV085OriginalPrimeFiniteGate
import WeilDefect.MasterReturn.CarrierTranslation
import Mathlib.Analysis.InnerProductSpace.Continuous

/-! RETURN109-F3: literal prime and cross-pole continuity on the
actual IV logarithmic carrier. No smoothness or graph membership premise.
The finite prime cutoff is exact, not a truncation error. -/
namespace WeilDefect.MasterReturn.LogCarrierPrimePole

open MeasureTheory InnerProductSpace
open scoped InnerProduct ComplexConjugate BigOperators
open CarrierTranslation
noncomputable section

def physical (a : ℝ) (f : NeutralLogHilbertCarrier a) : RealComplexL2 :=
  neutralLogPhysicalInclusion a f

theorem physical_continuous (a : ℝ) : Continuous (physical a) :=
  (neutralLogPhysicalInclusion a).continuous

theorem window_eq_ae (a : ℝ) (f : NeutralLogHilbertCarrier a) :
    neutralWindowRepresentative a (physical a f) =ᵐ[volume] physical a f := by
  filter_upwards [neutralLogHilbertCarrier_supported f] with x hx
  change neutralWindowRepresentative a (neutralLogPhysical f.val) x =
    neutralLogPhysical f.val x
  by_cases hm : x ∈ Set.Icc (-a) a
  · simp [neutralWindowRepresentative, hm]
  · simp [neutralWindowRepresentative, hm, hx hm]

theorem correlation_inner (a t : ℝ) (f g : NeutralLogHilbertCarrier a) :
    neutralWindowCorrelation a (physical a f) (physical a g) t =
      inner ℂ (shiftPhysical t (physical a f)) (physical a g) := by
  rw [L2.inner_def]
  unfold neutralWindowCorrelation
  apply integral_congr_ae
  filter_upwards [window_eq_ae a g,
    (shift_preserving t).quasiMeasurePreserving.ae (window_eq_ae a f),
    shiftPhysical_ae t (physical a f)] with x hg hf ht
  change neutralWindowRepresentative a (physical a f) (x-t) =
    physical a f (x-t) at hf
  change shiftPhysical t (physical a f) x = physical a f (x-t) at ht
  rw [hg, hf, ht]
  simp only [RCLike.inner_apply']
  ring

theorem shift_continuous (t : ℝ) :
    Continuous (shiftPhysical t : RealComplexL2 → RealComplexL2) :=
  (Lp.isometry_compMeasurePreserving (shift_preserving t)).continuous

theorem correlation_diagonal_continuous (a t : ℝ) :
    Continuous (fun f : NeutralLogHilbertCarrier a =>
      neutralWindowCorrelation a (physical a f) (physical a f) t) := by
  simp_rw [correlation_inner]
  exact ((shift_continuous t).comp (physical_continuous a)).inner
    (physical_continuous a)

theorem correlation_diagonal_norm_le (a t : ℝ) (f : NeutralLogHilbertCarrier a) :
    ‖neutralWindowCorrelation a (physical a f) (physical a f) t‖ ≤
      ‖physical a f‖ ^ 2 := by
  rw [correlation_inner]
  have h := norm_inner_le_norm (𝕜 := ℂ) (shiftPhysical t (physical a f)) (physical a f)
  have hs : ‖shiftPhysical t (physical a f)‖ = ‖physical a f‖ :=
    Lp.norm_compMeasurePreserving _ (shift_preserving t)
  rw [hs] at h
  simpa only [pow_two] using h

def primeFinite (a : ℝ) (N : ℕ) (f : NeutralLogHilbertCarrier a) : ℂ :=
  ∑ n ∈ Finset.range N, neutralIV085PrimeSummand
    (neutralWindowCorrelation a (physical a f) (physical a f)) n

theorem primeFinite_continuous (a : ℝ) (N : ℕ) :
    Continuous (primeFinite a N) := by
  unfold primeFinite
  apply continuous_finset_sum
  intro n hn
  unfold neutralIV085PrimeSummand
  exact continuous_const.mul
    ((correlation_diagonal_continuous a (Real.log n)).add
      (correlation_diagonal_continuous a (-Real.log n)))

theorem primeLiteral_eq_finite (a : ℝ) (N : ℕ)
    (hN : Real.exp (2*a) < (N : ℝ)) (f : NeutralLogHilbertCarrier a) :
    (∑' n : ℕ, neutralIV085PrimeSummand
      (neutralWindowCorrelation a (physical a f) (physical a f)) n) =
      primeFinite a N f := by
  apply neutralIV085PrimeSum_eq_finite _ _ N hN
  intro u hu
  apply abs_le.mpr
  by_contra h
  apply hu
  exact neutralWindowCorrelation_supported a _ _ u h

def poleDiagonal (a : ℝ) (f : NeutralLogHilbertCarrier a) : ℂ :=
  inner ℂ f (neutralLogPoleOperator a f)

theorem poleDiagonal_continuous (a : ℝ) : Continuous (poleDiagonal a) :=
  continuous_id.inner (neutralLogPoleOperator a).continuous

theorem poleLiteral_eq_diagonal (a : ℝ) (f : NeutralLogHilbertCarrier a) :
    Zeta23.paperFT (neutralWindowCorrelation a (physical a f) (physical a f))
      (Complex.I / 2) +
    Zeta23.paperFT (neutralWindowCorrelation a (physical a f) (physical a f))
      (-Complex.I / 2) = poleDiagonal a f := by
  have hp : (Complex.I / 2 : ℂ) = Complex.I * ((1/2 : ℝ) : ℂ) := by
    push_cast
    ring
  have hm : (-Complex.I / 2 : ℂ) = Complex.I * ((-(1/2) : ℝ) : ℂ) := by
    push_cast
    ring
  rw [hp, hm]
  exact neutralLogHilbertCorrelation_poleOperator a f f

theorem primePole_continuous (a : ℝ) (N : ℕ) :
    Continuous (fun f : NeutralLogHilbertCarrier a =>
      (poleDiagonal a f - primeFinite a N f).re) :=
  Complex.continuous_re.comp ((poleDiagonal_continuous a).sub
    (primeFinite_continuous a N))

end
end WeilDefect.MasterReturn.LogCarrierPrimePole

#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.physical_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.window_eq_ae
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.correlation_inner
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.shift_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.correlation_diagonal_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.correlation_diagonal_norm_le
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.primeFinite_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.primeLiteral_eq_finite
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.poleDiagonal_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.poleLiteral_eq_diagonal
#print axioms WeilDefect.MasterReturn.LogCarrierPrimePole.primePole_continuous
