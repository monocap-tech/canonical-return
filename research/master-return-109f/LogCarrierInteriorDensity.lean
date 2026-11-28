import WeilDefect.MasterReturn.LogCarrierMollifierLift
import Mathlib.Analysis.SpecificLimits.Normed

/-! F4c5c: assemble actual interior smooth form-core density.
All preceding staged declarations require native verification. -/
namespace WeilDefect.MasterReturn.LogCarrierInteriorDensity

open MeasureTheory Filter FourierTransform
open scoped Topology FourierTransform
open WeilDefect.MasterRecon
open LogCarrierCompressionLift LogCarrierMollifierLift
open LogCarrierDilationFourier LogCarrierInteriorMollifier
open LogCarrierFormDensityReceivers
noncomputable section

def interiorCore (a : ℝ) : Set (V08I.logCarrier a) :=
  {g | ∃ f : ℝ → ℂ, ContDiff ℝ ∞ f ∧ HasCompactSupport f ∧
    Function.support f ⊆ Set.Ioo (-a) a ∧
    (V08I.physicalInclusion a g : ℝ → ℂ) =ᵐ[volume] f}

theorem smoothLift_mem_core (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (φ : ContDiffBump (0 : ℝ)) (hr : s*a+φ.rOut < a)
    (g : V08I.logCarrier a) :
    smoothLift a s ha hs hs' φ hr g ∈ interiorCore a := by
  refine ⟨smoothWindow a s φ g, smoothWindow_smooth a s (by linarith) φ g,
    smoothWindow_compact a s ha (by linarith) φ g,
    smoothWindow_interior a s ha (by linarith) φ hr g, ?_⟩
  rw [smoothLift_physical]
  exact smoothL2_ae a s ha (by linarith) φ g

def radius (c : ℝ) (n : ℕ) : ℝ := c / (2*((n:ℝ)+1))

theorem radius_pos (c : ℝ) (hc : 0 < c) (n : ℕ) : 0 < radius c n := by
  unfold radius
  positivity

theorem radius_lt (c : ℝ) (hc : 0 < c) (n : ℕ) : radius c n < c := by
  unfold radius
  apply (div_lt_iff₀ (by positivity)).mpr
  nlinarith [Nat.cast_nonneg (α := ℝ) n,
    mul_nonneg hc.le (Nat.cast_nonneg (α := ℝ) n)]

theorem radius_tendsto (c : ℝ) : Tendsto (radius c) atTop (𝓝 0) := by
  have h : Tendsto (fun n : ℕ => (n:ℝ)+1) atTop atTop :=
    (tendsto_natCast_atTop_atTop (R := ℝ)).atTop_add_const_right 1
  have h' : Tendsto (fun n : ℕ => ((n:ℝ)+1)⁻¹) atTop (𝓝 0) :=
    tendsto_inv_atTop_zero.comp h
  have he : radius c = fun n : ℕ => (c/2)*((n:ℝ)+1)⁻¹ := by
    funext n
    unfold radius
    ring
  rw [he]
  simpa using (tendsto_const_nhds (x := c/2)).mul h'

def collarBump (c : ℝ) (hc : 0 < c) (n : ℕ) : ContDiffBump (0 : ℝ) :=
  ⟨radius c n/2, radius c n, by positivity [radius_pos c hc n],
    by linarith [radius_pos c hc n]⟩

theorem collarBump_rOut (c : ℝ) (hc : 0 < c) (n : ℕ) :
    (collarBump c hc n).rOut = radius c n := rfl

theorem smoothLift_tendsto (a s : ℝ) (ha : 0 ≤ a)
    (hs : (1/2 : ℝ) ≤ s) (hs' : s ≤ 1)
    (φ : ℕ → ContDiffBump (0 : ℝ))
    (hr : ∀ n, s*a+(φ n).rOut < a)
    (ht : Tendsto (fun n => (φ n).rOut) atTop (𝓝 0))
    (g : V08I.logCarrier a) :
    Tendsto (fun n => smoothLift a s ha hs hs' (φ n) (hr n) g)
      atTop (𝓝 (compression a s ha hs hs' g)) := by
  let m : ℕ → ℝ → ℂ := fun n => 𝓕 (bumpComplex (φ n))
  have hm (n : ℕ) : Measurable (m n) :=
    (VectorFourier.fourierIntegral_continuous (L := innerₗ ℝ)
      Real.continuous_fourierChar continuous_inner (bumpComplex_integrable (φ n))).measurable
  have hb (n : ℕ) (u : ℝ) : ‖m n u‖ ≤ 1 := bump_fourier_norm_le (φ n) u
  have hlim : ∀ᵐ u ∂volume, Tendsto (fun n => m n u) atTop (𝓝 1) :=
    Filter.Eventually.of_forall (bump_fourier_tendsto_one φ ht)
  have he := weighted_multiplier_error_integral_tendsto a
    (compression a s ha hs hs' g) m hm hb hlim
  simp only [compression_physical] at he
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩ := Metric.tendsto_atTop.mp he (ε^2) (by positivity)
  refine ⟨N, fun n hn => ?_⟩
  have h := hN n hn
  change dist (∫ u, V08I.logWeight u *
      ‖((𝓕 (bumpComplex (φ n))) u-1) *
        (𝓕 (LogCarrierCompressionL2.compressedL2 a s (by linarith) g) : V08I.FullL2) u‖^2) 0 <
      ε^2 at h
  rw [← smoothLift_error a s ha hs hs' (φ n) (hr n) g,
    Real.dist_eq, sub_zero, abs_of_nonneg (sq_nonneg _)] at h
  rw [dist_eq_norm]
  by_contra hbad
  have hh := le_of_not_gt hbad
  nlinarith [mul_nonneg (sub_nonneg.mpr hh)
    (show 0 ≤ ‖smoothLift a s ha hs hs' (φ n) (hr n) g -
      compression a s ha hs hs' g‖+ε by positivity)]

theorem interiorCore_dense (a : ℝ) (ha : 0 < a) :
    Dense (interiorCore a) := by
  intro g
  apply Metric.mem_closure_iff.mpr
  intro ε hε
  let s : ℕ → ℝ := fun n => 1-radius 1 n
  have hs (n : ℕ) : (1/2 : ℝ) ≤ s n := by
    dsimp [s, radius]
    have h : (1:ℝ)/(2*((n:ℝ)+1)) ≤ 1/2 := by
      apply (div_le_iff₀ (by positivity)).mpr
      linarith [Nat.cast_nonneg (α := ℝ) n]
    linarith
  have hslt (n : ℕ) : s n < 1 := by
    dsimp [s]
    linarith [radius_pos 1 (by norm_num) n]
  have hs' (n : ℕ) : s n ≤ 1 := (hslt n).le
  have ht : Tendsto s atTop (𝓝 1) := by
    simpa [s] using (tendsto_const_nhds (x := (1:ℝ))).sub (radius_tendsto 1)
  obtain ⟨N,hN⟩ := Metric.tendsto_atTop.mp
    (compression_tendsto a ha.le g s hs hs' ht) (ε/2) (by positivity)
  let t := s N
  let c := a-t*a
  have hc : 0 < c := by
    dsimp [c,t]
    nlinarith [mul_pos ha (sub_pos.mpr (hslt N))]
  let φ := collarBump c hc
  have hr (n : ℕ) : t*a+(φ n).rOut < a := by
    change t*a+radius c n < a
    linarith [radius_lt c hc n]
  have hφ : Tendsto (fun n => (φ n).rOut) atTop (𝓝 0) := radius_tendsto c
  obtain ⟨M,hM⟩ := Metric.tendsto_atTop.mp
    (smoothLift_tendsto a t ha.le (hs N) (hs' N) φ hr hφ g)
    (ε/2) (by positivity)
  let q := smoothLift a t ha.le (hs N) (hs' N) (φ M) (hr M) g
  refine ⟨q, smoothLift_mem_core a t ha.le (hs N) (hs' N) (φ M) (hr M) g, ?_⟩
  have h1 := hN N le_rfl
  have h2 := hM M le_rfl
  rw [dist_comm g q]
  exact (dist_triangle q (compression a t ha.le (hs N) (hs' N) g) g).trans_lt
    (by linarith)

end
end WeilDefect.MasterReturn.LogCarrierInteriorDensity

#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.smoothLift_mem_core
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.radius_pos
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.radius_lt
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.radius_tendsto
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.collarBump_rOut
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.smoothLift_tendsto
#print axioms WeilDefect.MasterReturn.LogCarrierInteriorDensity.interiorCore_dense
