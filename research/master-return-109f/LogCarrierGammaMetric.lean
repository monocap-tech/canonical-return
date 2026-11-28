import WeilDefect.MasterRecon.V08I
import WeilDefect.Arithmetic.ActualZetaIV071GammaFinitePrelude
import Mathlib.Tactic

/-! RETURN109-F2: scalar comparison in the exact receiving metric.
The global logarithmic gain estimate and measurability remain explicit
premises until separately proved. No replacement of the actual symbol. -/
namespace WeilDefect.MasterReturn.LogCarrierGammaMetric

open MeasureTheory
noncomputable section

theorem angular_log_gap (r : ℝ) :
    0 ≤ Real.log (2 + |r|) - (1 / 2 : ℝ) * Real.log (1 + r^2) ∧
    Real.log (2 + |r|) - (1 / 2 : ℝ) * Real.log (1 + r^2) ≤
      (1 / 2 : ℝ) * Real.log 5 := by
  have hp : 0 < 2 + |r| := by positivity
  have hq : 0 < 1 + r^2 := by positivity
  have hs : |r|^2 = r^2 := sq_abs r
  have hlo : 1 + r^2 ≤ (2 + |r|)^2 := by
    nlinarith [hs, abs_nonneg r]
  have hhi : (2 + |r|)^2 ≤ 5 * (1 + r^2) := by
    nlinarith [hs, sq_nonneg (2 * |r| - 1)]
  have hl := Real.log_le_log hq hlo
  have hsq : 0 < (2 + |r|)^2 := by positivity
  have hu := Real.log_le_log hsq hhi
  rw [Real.log_pow] at hl hu
  rw [Real.log_mul (by norm_num : (5 : ℝ) ≠ 0) hq.ne'] at hu
  norm_num at hl hu
  constructor <;> linarith

def LogGainBound : Prop :=
  ∀ r : ℝ,
    0 ≤ Zeta23.EF.gammaBracket r - Zeta23.EF.gammaBracket 0 -
      (1 / 2 : ℝ) * Real.log (1 + r^2) ∧
    Zeta23.EF.gammaBracket r - Zeta23.EF.gammaBracket 0 -
      (1 / 2 : ℝ) * Real.log (1 + r^2) ≤ 4 + Real.log 2

theorem receiving_remainder_bounds (h : LogGainBound) (u : ℝ) :
    Zeta23.EF.gammaBracket 0 - (1 / 2 : ℝ) * Real.log 5 ≤
      Zeta23.EF.gammaBracket (2 * Real.pi * u) - V08I.logWeight u ∧
    Zeta23.EF.gammaBracket (2 * Real.pi * u) - V08I.logWeight u ≤
      Zeta23.EF.gammaBracket 0 + 4 + Real.log 2 := by
  obtain ⟨hl, hu⟩ := h (2 * Real.pi * u)
  obtain ⟨gl, gu⟩ := angular_log_gap (2 * Real.pi * u)
  unfold V08I.logWeight
  constructor <;> linarith

def remainderBound : ℝ :=
  |Zeta23.EF.gammaBracket 0| + 4 + Real.log 2 + (1 / 2 : ℝ) * Real.log 5

theorem remainderBound_nonneg : 0 ≤ remainderBound := by
  unfold remainderBound
  have h2 := Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 2)
  have h5 := Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 5)
  linarith [abs_nonneg (Zeta23.EF.gammaBracket 0)]

theorem receiving_remainder_abs (h : LogGainBound) (u : ℝ) :
    |Zeta23.EF.gammaBracket (2 * Real.pi * u) - V08I.logWeight u| ≤
      remainderBound := by
  obtain ⟨hl, hu⟩ := receiving_remainder_bounds h u
  have h2 := Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 2)
  have h5 := Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 5)
  unfold remainderBound
  apply abs_le.mpr
  constructor <;> linarith [le_abs_self (Zeta23.EF.gammaBracket 0),
    neg_abs_le (Zeta23.EF.gammaBracket 0)]

def gammaRemainder (u : ℝ) : ℂ :=
  ((Zeta23.EF.gammaBracket (2 * Real.pi * u) - V08I.logWeight u : ℝ) : ℂ)

theorem gammaRemainder_norm (h : LogGainBound) (u : ℝ) :
    ‖gammaRemainder u‖ ≤ remainderBound := by
  rw [gammaRemainder, Complex.norm_real, Real.norm_eq_abs]
  exact receiving_remainder_abs h u

def gammaRemainderOperator (h : LogGainBound)
    (hm : Measurable gammaRemainder) : V08I.FullL2 →L[ℂ] V08I.FullL2 :=
  V08I.multiplier gammaRemainder hm remainderBound (gammaRemainder_norm h)

theorem gammaRemainderOperator_ae (h : LogGainBound)
    (hm : Measurable gammaRemainder) (f : V08I.FullL2) :
    (gammaRemainderOperator h hm f : ℝ → ℂ) =ᵐ[volume]
      fun u => gammaRemainder u * f u :=
  V08I.multiplier_ae gammaRemainder hm remainderBound (gammaRemainder_norm h) f

end
end WeilDefect.MasterReturn.LogCarrierGammaMetric

#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.angular_log_gap
#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.receiving_remainder_bounds
#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.remainderBound_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.receiving_remainder_abs
#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.gammaRemainder_norm
#print axioms WeilDefect.MasterReturn.LogCarrierGammaMetric.gammaRemainderOperator_ae
