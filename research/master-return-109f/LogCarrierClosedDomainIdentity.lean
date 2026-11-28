import WeilDefect.MasterReturn.LogCarrierBoundedAnalysis
import WeilDefect.MasterReturn.LogCarrierClosedDomainCore
import WeilDefect.MasterReturn.LogCarrierPrimePole
import WeilDefect.MasterReturn.LogCarrierArchContinuity
import WeilDefect.MasterReturn.LogCarrierInteriorDensity
import Mathlib.Topology.Algebra.Module.Equiv

/-! RETURN109-F5: actual closed-domain arithmetic/graph equality.
Candidate only: imports include E and the independent F proof banks.
Arithmetic side is the literal literature RHS, not defined via signedGram.
No positivity, rigidity or additional smoothness hypothesis. -/
namespace WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity

open MeasureTheory InnerProductSpace
open scoped InnerProduct
open WeilDefect.MasterRecon
open LogCarrierContinuous LogCarrierBoundedAnalysis LogCarrierClosedDomainCore
open LogCarrierPrimePole LogCarrierArchContinuity LogCarrierInteriorDensity
open LogCarrierGlobalEnergy
noncomputable section
set_option maxHeartbeats 3000000
set_option maxRecDepth 4096

def arithmeticForm (a : ℝ) (g : NeutralLogHilbertCarrier a) : ℝ :=
  (Zeta23.EF.literatureRHS
    (neutralWindowCorrelation a (neutralLogPhysical g.val)
      (neutralLogPhysical g.val))).re

def signedForm (a : ℝ)
    (T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index)
    (g : NeutralLogHilbertCarrier a) : ℝ :=
  (inner ℂ (V10B6c.signedGram T (carrierContinuousEquiv a g))
    (carrierContinuousEquiv a g)).re

def ivInteriorCore (a : ℝ) : Set (NeutralLogHilbertCarrier a) :=
  (carrierContinuousEquiv a) ⁻¹' interiorCore a

theorem ivInteriorCore_dense (a : ℝ) (ha : 0 < a) :
    Dense (ivInteriorCore a) :=
  (interiorCore_dense a ha).preimage (carrierContinuousEquiv a).toHomeomorph.isOpenMap

theorem ivCore_window_smooth (a : ℝ) (g : NeutralLogHilbertCarrier a)
    (hg : g ∈ ivInteriorCore a) :
    ∃ ψ : ℝ → ℂ, ContDiff ℝ 2 ψ ∧ HasCompactSupport ψ ∧
      neutralWindowRepresentative a (neutralLogPhysical g.val) =ᵐ[volume] ψ := by
  rcases hg with ⟨ψ,hψ,hc,hs,hp⟩
  have he : (V08I.physicalInclusion a (carrierContinuousEquiv a g) : ℝ → ℂ) =
      (neutralLogPhysical g.val : ℝ → ℂ) := by rw [carrierContinuousEquiv_physical]
  rw [he] at hp
  refine ⟨ψ,hψ.of_le (by norm_num),hc,?_⟩
  exact (window_eq_ae a g).trans hp

theorem arithmeticForm_decompose (a : ℝ) (N : ℕ)
    (hN : Real.exp (2*a) < (N : ℝ)) (g : NeutralLogHilbertCarrier a) :
    arithmeticForm a g =
      (poleDiagonal a g - primeFinite a N g +
        literalArch a (carrierContinuousEquiv a g)).re := by
  unfold arithmeticForm Zeta23.EF.literatureRHS
  have hphys : V08I.physicalInclusion a (carrierContinuousEquiv a g) = physical a g :=
    carrierContinuousEquiv_physical a g
  unfold literalArch
  rw [hphys]
  change (Zeta23.paperFT (neutralWindowCorrelation a (physical a g) (physical a g))
    (Complex.I/2) +
    Zeta23.paperFT (neutralWindowCorrelation a (physical a g) (physical a g))
    (-Complex.I/2) -
    (∑' n : ℕ, neutralIV085PrimeSummand
      (neutralWindowCorrelation a (physical a g) (physical a g)) n) +
    (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ,
      Zeta23.paperFT (neutralWindowCorrelation a (physical a g) (physical a g)) (r : ℂ) *
      (Zeta23.EF.gammaBracket r : ℂ)).re = _
  rw [poleLiteral_eq_diagonal, primeLiteral_eq_finite a N hN]

theorem arithmeticForm_continuous (a : ℝ) : Continuous (arithmeticForm a) := by
  obtain ⟨N,hN⟩ := exists_nat_gt (Real.exp (2*a))
  have he : arithmeticForm a = fun g =>
      (poleDiagonal a g - primeFinite a N g +
        literalArch a (carrierContinuousEquiv a g)).re := by
    funext g
    exact arithmeticForm_decompose a N hN g
  rw [he]
  exact Complex.continuous_re.comp
    (((poleDiagonal_continuous a).sub (primeFinite_continuous a N)).add
      ((literalArch_continuous a).comp (carrierContinuousEquiv a).continuous))

theorem signedForm_continuous (a : ℝ)
    (T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index) :
    Continuous (signedForm a T) :=
  Complex.continuous_re.comp
    (((V10B6c.signedGram T).continuous.comp (carrierContinuousEquiv a).continuous).inner
      (𝕜 := ℂ) (carrierContinuousEquiv a).continuous)

theorem graph_core_identity (a : ℝ) (h : NeutralActualZetaSourceDomain a)
    (hc : neutralActualZetaSourcePhysical a h ∈ ivInteriorCore a) :
    neutralActualZetaSourceQuadratic a h =
      arithmeticForm a (neutralActualZetaSourcePhysical a h) := by
  obtain ⟨ψ,hψ,hcomp,hae⟩ := ivCore_window_smooth a _ hc
  exact smooth_graph_literature a h ψ hψ hcomp hae

theorem graph_signed_of_literal (a : ℝ)
    (T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index)
    (hT : V10B6d.LiteralCoordinates a T) (h : NeutralActualZetaSourceDomain a) :
    neutralActualZetaSourceQuadratic a h =
      signedForm a T (neutralActualZetaSourcePhysical a h) := by
  unfold signedForm
  rw [graph_signed_energy, V10B6c.signedGram_diagonal,
    ← V10B6d.signed_channel_tsum (T (carrierContinuousEquiv a
      (neutralActualZetaSourcePhysical a h)))]
  apply tsum_congr
  intro i
  rw [hT (carrierContinuousEquiv a (neutralActualZetaSourcePhysical a h)) i]
  rfl

theorem actual_closed_domain_identity (a : ℝ) (ha : 0 < a) :
    ∃ T : V08I.logCarrier a →L[ℂ] V07C.Coeff V10A.Index,
      V10B6d.LiteralCoordinates a T ∧
      (∀ g : NeutralLogHilbertCarrier a, arithmeticForm a g = signedForm a T g) ∧
      (∀ g : NeutralLogHilbertCarrier a, ∃ h : NeutralActualZetaSourceDomain a,
        neutralActualZetaSourcePhysical a h = g ∧
        graphCoeff a h = T (carrierContinuousEquiv a g) ∧
        neutralActualZetaSourceQuadratic a h = arithmeticForm a g) ∧
      (∀ h : NeutralActualZetaSourceDomain a,
        neutralActualZetaSourceQuadratic a h =
          arithmeticForm a (neutralActualZetaSourcePhysical a h)) := by
  obtain ⟨T,hT,hgraph⟩ := actual_graph_receiver a ha.le
  have he : ∀ g, arithmeticForm a g = signedForm a T g := by
    apply continuous_extension_unique a (ivInteriorCore a) (ivInteriorCore_dense a ha)
      _ _ (arithmeticForm_continuous a) (signedForm_continuous a T)
    intro g hg
    obtain ⟨h,hp,hcoeff,hq⟩ := hgraph g
    have hc : neutralActualZetaSourcePhysical a h ∈ ivInteriorCore a := hp.symm ▸ hg
    have hs := graph_core_identity a h hc
    rw [hp] at hs
    exact hs.symm.trans hq
  refine ⟨T,hT,he,?_,?_⟩
  · intro g
    obtain ⟨h,hp,hcoeff,hq⟩ := hgraph g
    exact ⟨h,hp,hcoeff,hq.trans (he g).symm⟩
  · intro h
    exact (graph_signed_of_literal a T hT h).trans
      (he (neutralActualZetaSourcePhysical a h)).symm

end
end WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity

#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.ivInteriorCore_dense
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.ivCore_window_smooth
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.arithmeticForm_decompose
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.arithmeticForm_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.signedForm_continuous
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.graph_core_identity
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.graph_signed_of_literal
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainIdentity.actual_closed_domain_identity
