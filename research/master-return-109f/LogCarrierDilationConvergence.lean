import WeilDefect.MasterReturn.LogCarrierDilationEnergy
import Mathlib.MeasureTheory.Integral.DominatedConvergence

/-! F4c4c: strong weighted dilation convergence by compact/tail
splitting. General continuous bounded ordinary Fourier coordinate
with integrable logarithmic energy; no weighted-space density axiom. -/
namespace WeilDefect.MasterReturn.LogCarrierDilationConvergence

open MeasureTheory Filter
open scoped Topology FourierTransform
open WeilDefect.MasterRecon
open LogCarrierDilationEnergy LogCarrierCompressionL2
noncomputable section

def err (s : ℝ) (F : ℝ → ℂ) (u : ℝ) : ℝ :=
  energy (fun v => (s : ℂ)*F (s*v)-F v) u

theorem err_nonneg (s : ℝ) (F : ℝ → ℂ) (u : ℝ) : 0 ≤ err s F u :=
  energy_nonneg _ _

theorem err_le (s : ℝ) (F : ℝ → ℂ) (u : ℝ) :
    err s F u ≤ 2*dilatedEnergy s F u + 2*energy F u := by
  let z := (s : ℂ)*F (s*u)
  have h := norm_sub_le z (F u)
  have hh : ‖z-F u‖^2 ≤ 2*‖z‖^2+2*‖F u‖^2 := by
    nlinarith [sq_nonneg (‖z‖-‖F u‖),
      mul_nonneg (sub_nonneg.mpr h) (show 0 ≤ ‖z‖+‖F u‖+‖z-F u‖ by positivity)]
  unfold err energy dilatedEnergy
  change V08I.logWeight u * ‖z-F u‖^2 ≤ _
  have hh' := mul_le_mul_of_nonneg_left hh (V08I.logWeight_pos u).le
  dsimp [z] at hh'
  nlinarith only [hh']

theorem err_integrable (s : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (F : ℝ → ℂ) (hF : Continuous F) (hE : Integrable (energy F)) :
    Integrable (err s F) := by
  have h := ((dilatedEnergy_integrable s hs hs' F hF hE).const_mul 2).add
    (hE.const_mul 2)
  apply h.mono' (show Continuous (err s F) from by unfold err energy; fun_prop).aestronglyMeasurable
  filter_upwards [] with u
  rw [Real.norm_eq_abs, abs_of_nonneg (err_nonneg s F u)]
  exact err_le s F u

theorem compact_error_tendsto (F : ℝ → ℂ) (hF : Continuous F)
    (B : ℝ) (hB : 0 ≤ B) (hb : ∀ u, ‖F u‖ ≤ B)
    (s : ℕ → ℝ) (hs : ∀ n, (1/2 : ℝ) ≤ s n) (hs' : ∀ n, s n ≤ 1)
    (ht : Tendsto s atTop (𝓝 1)) (R : ℝ) :
    Tendsto (fun n => ∫ u in Set.Icc (-R) R, err (s n) F u) atTop (𝓝 0) := by
  let μ := volume.restrict (Set.Icc (-R) R)
  have hd : Integrable (fun u => 4*B^2*V08I.logWeight u) μ :=
    (show Continuous (fun u => 4*B^2*V08I.logWeight u) from by fun_prop).integrableOn_Icc
  have h := tendsto_integral_of_dominated_convergence (μ := μ)
    (F := fun n u => err (s n) F u) (f := fun _ => (0 : ℝ))
    (fun u => 4*B^2*V08I.logWeight u)
    (fun n => (show Continuous (err (s n) F) from by unfold err energy; fun_prop).aestronglyMeasurable)
    hd
    (fun n => Filter.Eventually.of_forall fun u => by
      have hs0 : 0 ≤ s n := by linarith [hs n]
      have hz : ‖(s n : ℂ)*F (s n*u)‖ ≤ B := by
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hs0]
        calc
          _ ≤ s n * B := mul_le_mul_of_nonneg_left (hb _) hs0
          _ ≤ B := by nlinarith [mul_nonneg hB (sub_nonneg.mpr (hs' n))]
      have hn := norm_sub_le ((s n : ℂ)*F (s n*u)) (F u)
      have h2 : ‖(s n : ℂ)*F (s n*u)-F u‖ ≤ 2*B := by linarith [hb u]
      have hsq : ‖(s n : ℂ)*F (s n*u)-F u‖^2 ≤ 4*B^2 := by
        nlinarith [norm_nonneg ((s n : ℂ)*F (s n*u)-F u),
          mul_nonneg (sub_nonneg.mpr h2)
            (show 0 ≤ 2*B+‖(s n : ℂ)*F (s n*u)-F u‖ by positivity)]
      rw [Real.norm_eq_abs, abs_of_nonneg (err_nonneg (s n) F u)]
      unfold err energy
      have hh := mul_le_mul_of_nonneg_left hsq (V08I.logWeight_pos u).le
      nlinarith only [hh])
    (Filter.Eventually.of_forall fun u => by
      have hu : Tendsto (fun t : ℝ => err t F u) (𝓝 1) (𝓝 (err 1 F u)) :=
        (show Continuous (fun t : ℝ => err t F u) from by unfold err energy; fun_prop).continuousAt
      have he : err 1 F u = 0 := by simp [err, energy]
      rw [he] at hu
      exact hu.comp ht)
  simpa [μ] using h

theorem tail_energy_tendsto (F : ℝ → ℂ) (hE : Integrable (energy F)) :
    Tendsto (fun N : ℕ => ∫ u in tail (N : ℝ), energy F u) atTop (𝓝 0) := by
  have hanti : Antitone (fun N : ℕ => tail (N : ℝ)) := by
    intro i j hij u hu
    exact lt_of_le_of_lt (by exact_mod_cast hij) hu
  have h := tendsto_setIntegral_of_antitone (μ := volume) (f := energy F)
    (fun N : ℕ => tail_measurable (N : ℝ)) hanti ⟨0,hE.integrableOn⟩
  have he : (⋂ N : ℕ, tail (N : ℝ)) = ∅ := by
    apply Set.eq_empty_iff_forall_notMem.mpr
    intro u hu
    obtain ⟨N,hN⟩ := exists_nat_gt |u|
    exact (not_lt_of_ge hN.le) (Set.mem_iInter.mp hu N)
  rw [he, integral_empty] at h
  exact h

theorem err_tail_budget (s R : ℝ) (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (hR : 0 ≤ R) (F : ℝ → ℂ) (hF : Continuous F) (hE : Integrable (energy F)) :
    (∫ u in tail R, err s F u) ≤ 6 * ∫ u in tail (R/2), energy F u := by
  have hd := dilatedEnergy_integrable s hs hs' F hF hE
  have hh := integral_mono_ae
    ((err_integrable s hs hs' F hF hE).integrableOn : IntegrableOn (err s F) (tail R))
    (((hd.const_mul 2).add (hE.const_mul 2)).integrableOn :
      IntegrableOn (fun u => 2*dilatedEnergy s F u+2*energy F u) (tail R))
    (Filter.Eventually.of_forall (err_le s F))
  rw [integral_add (hd.const_mul 2).integrableOn (hE.const_mul 2).integrableOn,
    integral_const_mul, integral_const_mul] at hh
  have hb := dilatedEnergy_tail s R hs hs' hR F hE
  have hm : (∫ u in tail R, energy F u) ≤ ∫ u in tail (R/2), energy F u := by
    apply integral_mono_measure
    · apply Measure.restrict_mono
      · intro u hu
        change R < |u| at hu
        change R/2 < |u|
        linarith
      · rfl
    · exact Filter.Eventually.of_forall (energy_nonneg F)
    · exact hE.integrableOn
  have hpos : 0 ≤ ∫ u in tail (R/2), energy F u :=
    integral_nonneg (energy_nonneg F)
  nlinarith [mul_nonneg (sub_nonneg.mpr hs') hpos]

theorem strong_weighted_dilation (F : ℝ → ℂ) (hF : Continuous F)
    (B : ℝ) (hB : 0 ≤ B) (hb : ∀ u, ‖F u‖ ≤ B) (hE : Integrable (energy F))
    (s : ℕ → ℝ) (hs : ∀ n, (1/2 : ℝ) ≤ s n) (hs' : ∀ n, s n ≤ 1)
    (ht : Tendsto s atTop (𝓝 1)) :
    Tendsto (fun n => ∫ u, err (s n) F u) atTop (𝓝 0) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩ := Metric.tendsto_atTop.mp (tail_energy_tendsto F hE)
    (ε/12) (by positivity)
  have htail : (∫ u in tail (N : ℝ), energy F u) < ε/12 := by
    have h := hN N le_rfl
    rw [Real.dist_eq, sub_zero, abs_of_nonneg (integral_nonneg (energy_nonneg F))] at h
    exact h
  let R : ℝ := 2*(N : ℝ)
  have hR : 0 ≤ R := by dsimp [R]; positivity
  obtain ⟨M,hM⟩ := Metric.tendsto_atTop.mp
    (compact_error_tendsto F hF B hB hb s hs hs' ht R) (ε/2) (by positivity)
  refine ⟨M, fun n hn => ?_⟩
  have hc := hM n hn
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (integral_nonneg (err_nonneg (s n) F))] at hc
  have hbd := err_tail_budget (s n) R (hs n) (hs' n) hR F hF hE
  have hhalf : R/2 = (N : ℝ) := by dsimp [R]; ring
  rw [hhalf] at hbd
  have he : (Set.Icc (-R) R)ᶜ = tail R := by
    ext u
    simp only [Set.mem_compl_iff, Set.mem_Icc, tail, Set.mem_setOf_eq, ← abs_le, not_le]
  have hsum := integral_add_compl (s := Set.Icc (-R) R) measurableSet_Icc
    (err_integrable (s n) (hs n) (hs' n) F hF hE)
  rw [he] at hsum
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (integral_nonneg (err_nonneg (s n) F))]
  linarith

end
end WeilDefect.MasterReturn.LogCarrierDilationConvergence

#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.err_nonneg
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.err_le
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.err_integrable
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.compact_error_tendsto
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.tail_energy_tendsto
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.err_tail_budget
#print axioms WeilDefect.MasterReturn.LogCarrierDilationConvergence.strong_weighted_dilation
