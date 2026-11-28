import WeilDefect.MasterRecon.V08I
import WeilDefect.External.Zeta23.Theorems.Thm_Zeta23_WeilEF_digamma_growth_strip
import Mathlib.Tactic

/-! F2b: derive the weaker logarithmic envelope needed for F
from the actual existing strip-growth supplier. It does not imply
the coefficient-one bounded remainder needed for G. -/
namespace WeilDefect.MasterReturn.LogCarrierGammaGrowth

noncomputable section

theorem gamma_logWeight_bound :
    ∃ K : ℝ, 0 ≤ K ∧ ∀ u : ℝ,
      |Zeta23.EF.gammaBracket (2 * Real.pi * u)| ≤ K * V08I.logWeight u := by
  obtain ⟨C, hC, h⟩ := Zeta23.WeilEF.digamma_growth_strip
  have hlog : 0 < Real.log 2 := Real.log_pos (by norm_num)
  refine ⟨C + |Real.log Real.pi| / Real.log 2, by positivity, ?_⟩
  intro u
  let r : ℝ := 2 * Real.pi * u
  let z : ℂ := 1 / 4 + Complex.I * (r : ℂ) / 2
  have hre : z.re = 1 / 4 := by simp [z]
  have him : z.im = r / 2 := by simp [z]
  have hd := h z (by rw [hre]) (by rw [hre]; norm_num)
  have hs : |r / 2| ≤ |r| := by
    rw [abs_div, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 2)]
    linarith [abs_nonneg r]
  have hw : Real.log (2 + |z.im|) ≤ V08I.logWeight u := by
    rw [him]
    change Real.log (2 + |r / 2|) ≤ Real.log (2 + |r|)
    exact Real.log_le_log (by positivity) (by linarith)
  have hconst : |Real.log Real.pi| ≤
      (|Real.log Real.pi| / Real.log 2) * V08I.logWeight u := by
    calc
      _ = (|Real.log Real.pi| / Real.log 2) * Real.log 2 := by
        field_simp [hlog.ne']
      _ ≤ _ := mul_le_mul_of_nonneg_left (V08I.logWeight_lower u) (by positivity)
  have hg : Zeta23.EF.gammaBracket r = (Complex.digamma z).re - Real.log Real.pi := rfl
  change |Zeta23.EF.gammaBracket r| ≤ _
  rw [hg]
  calc
    _ ≤ |(Complex.digamma z).re| + |Real.log Real.pi| := abs_sub _ _
    _ ≤ ‖Complex.digamma z‖ + |Real.log Real.pi| :=
      add_le_add_right (Complex.abs_re_le_norm _) _
    _ ≤ C * V08I.logWeight u + (|Real.log Real.pi| / Real.log 2) *
        V08I.logWeight u :=
      add_le_add (hd.trans (mul_le_mul_of_nonneg_left hw hC.le)) hconst
    _ = _ := by ring

end
end WeilDefect.MasterReturn.LogCarrierGammaGrowth

#print axioms WeilDefect.MasterReturn.LogCarrierGammaGrowth.gamma_logWeight_bound
