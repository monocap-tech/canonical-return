import WeilDefect.Arithmetic.ActualZetaFullGraphRawDivisorQuadratic
import WeilDefect.Arithmetic.ActualZetaSmoothCorrelationRegularity

/-! RETURN109-F1 staging: exact a.e.-smooth original graph core and
uniqueness of continuous extension on the actual logarithmic carrier.
This module does not assert density, the Gamma multiplier estimate,
universal graph coverage, or an unconditional closed-domain identity. -/
namespace WeilDefect.MasterReturn.LogCarrierClosedDomainCore

open MeasureTheory
noncomputable section

theorem smooth_graph_literature (a : ℝ)
    (h : NeutralActualZetaSourceDomain a) (ψ : ℝ → ℂ)
    (hψ : ContDiff ℝ 2 ψ) (hc : HasCompactSupport ψ)
    (hae : neutralWindowRepresentative a
      (neutralLogPhysical (neutralActualZetaSourcePhysical a h).val)
        =ᵐ[volume] ψ) :
    neutralActualZetaSourceQuadratic a h =
      (Zeta23.EF.literatureRHS
        (neutralWindowCorrelation a
          (neutralLogPhysical (neutralActualZetaSourcePhysical a h).val)
          (neutralLogPhysical (neutralActualZetaSourcePhysical a h).val))).re := by
  apply neutralActualZetaSourceQuadratic_fullGraph_literature_of_contDiff
  exact neutralWindowCorrelation_contDiff_of_ae_smoothSecond a _ _ ψ hψ hc hae

theorem continuous_extension_unique (a : ℝ)
    (S : Set (NeutralLogHilbertCarrier a)) (hS : Dense S)
    (Q R : NeutralLogHilbertCarrier a → ℝ)
    (hQ : Continuous Q) (hR : Continuous R)
    (hcore : ∀ g ∈ S, Q g = R g) :
    ∀ g, Q g = R g := by
  have hs : S ⊆ {g | Q g = R g} := hcore
  have hc : closure S ⊆ {g | Q g = R g} :=
    closure_minimal hs (isClosed_eq hQ hR)
  intro g
  apply hc
  rw [hS.closure_eq]
  exact Set.mem_univ g

end
end WeilDefect.MasterReturn.LogCarrierClosedDomainCore

#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainCore.smooth_graph_literature
#print axioms WeilDefect.MasterReturn.LogCarrierClosedDomainCore.continuous_extension_unique
