import WeilDefect.MasterReturn.LogCarrierContinuous

/-! RETURN108-Y: actual isometric aperture inclusions and naturality of the
continuous carrier equivalence. No assertion about translation or signed forms. -/
namespace WeilDefect.MasterReturn.LogCarrierNatural
open MeasureTheory Set
open WeilDefect.MasterRecon
open LogCarrierContinuous

theorem iv_mem_of_le (a b : ℝ) (hab : a ≤ b) (g : NeutralLogHilbertCarrier a) :
    g.val ∈ NeutralLogHilbertCarrier b := by
  change neutralOutsideRestriction b (neutralLogPhysical g.val) = 0
  apply (neutralOutsideRestriction_zero_iff b _).mpr
  have hsub : Set.Icc (-a) a ⊆ Set.Icc (-b) b :=
    Set.Icc_subset_Icc (neg_le_neg hab) hab
  filter_upwards [neutralLogHilbertCarrier_supported g] with x hx
  intro hout
  exact hx (fun hin => hout (hsub hin))

theorem recon_mem_of_le (a b : ℝ) (hab : a ≤ b) (g : V08I.logCarrier a) :
    g.val ∈ V08I.logCarrier b := by
  change V08I.supportConstraint b (g : V08I.FullL2) = 0
  apply Lp.ext
  have hm := V08I.multiplier_ae (V08I.outsideMask b) (V08I.outsideMask_measurable b) 1
    (V08I.outsideMask_bound b) (V08I.physicalMap (g : V08I.FullL2))
  change (V08I.supportConstraint b (g : V08I.FullL2) : ℝ → ℂ) =ᵐ[volume]
    fun x => V08I.outsideMask b x * V08I.physicalMap (g : V08I.FullL2) x at hm
  have hsub : Set.Icc (-a) a ⊆ Set.Icc (-b) b :=
    Set.Icc_subset_Icc (neg_le_neg hab) hab
  filter_upwards [hm, V08I.physicalInclusion_support a g,
    Lp.coeFn_zero ℂ 2 (volume : Measure ℝ)] with x hx hs h0
  change V08I.supportConstraint b (g : V08I.FullL2) x = (0 : V08I.FullL2) x
  change (0 : V08I.FullL2) x = 0 at h0
  rw [hx, h0, V08I.outsideMask]
  by_cases h : x ∈ Set.Icc (-b) b
  · have hn : x ∉ (Set.Icc (-b) b)ᶜ := by simpa using h
    rw [Set.indicator_of_notMem hn, zero_mul]
  · have hz : V08I.physicalMap (g : V08I.FullL2) x = 0 :=
      hs (fun hin => h (hsub hin))
    rw [Set.indicator_of_mem (show x ∈ (Set.Icc (-b) b)ᶜ from h), hz, mul_zero]

noncomputable def ivInclLinear (a b : ℝ) (hab : a ≤ b) :
    NeutralLogHilbertCarrier a →ₗ[ℂ] NeutralLogHilbertCarrier b where
  toFun g := ⟨g.val, iv_mem_of_le a b hab g⟩
  map_add' _ _ := Subtype.ext rfl
  map_smul' _ _ := Subtype.ext rfl

noncomputable def reconInclLinear (a b : ℝ) (hab : a ≤ b) :
    V08I.logCarrier a →ₗ[ℂ] V08I.logCarrier b where
  toFun g := ⟨g.val, recon_mem_of_le a b hab g⟩
  map_add' _ _ := Subtype.ext rfl
  map_smul' _ _ := Subtype.ext rfl

noncomputable def ivIncl (a b : ℝ) (hab : a ≤ b) :
    NeutralLogHilbertCarrier a →L[ℂ] NeutralLogHilbertCarrier b :=
  (ivInclLinear a b hab).mkContinuous 1 (by
    intro g
    change ‖g‖ ≤ 1 * ‖g‖
    simp)

noncomputable def reconIncl (a b : ℝ) (hab : a ≤ b) :
    V08I.logCarrier a →L[ℂ] V08I.logCarrier b :=
  (reconInclLinear a b hab).mkContinuous 1 (by
    intro g
    change ‖g‖ ≤ 1 * ‖g‖
    simp)

theorem ivIncl_norm (a b : ℝ) (hab : a ≤ b) (g : NeutralLogHilbertCarrier a) :
    ‖ivIncl a b hab g‖ = ‖g‖ := rfl

theorem reconIncl_norm (a b : ℝ) (hab : a ≤ b) (g : V08I.logCarrier a) :
    ‖reconIncl a b hab g‖ = ‖g‖ := rfl

theorem ivIncl_physical (a b : ℝ) (hab : a ≤ b) (g : NeutralLogHilbertCarrier a) :
    neutralLogPhysical (ivIncl a b hab g).val = neutralLogPhysical g.val := rfl

theorem reconIncl_physical (a b : ℝ) (hab : a ≤ b) (g : V08I.logCarrier a) :
    V08I.physicalInclusion b (reconIncl a b hab g) =
      V08I.physicalInclusion a g := rfl

theorem carrier_natural (a b : ℝ) (hab : a ≤ b) (g : NeutralLogHilbertCarrier a) :
    carrierContinuousEquiv b (ivIncl a b hab g) =
      reconIncl a b hab (carrierContinuousEquiv a g) := by
  apply V08J.physicalInclusion_injective b
  rw [carrierContinuousEquiv_physical, ivIncl_physical, reconIncl_physical,
    carrierContinuousEquiv_physical]

theorem ivIncl_comp (a b c : ℝ) (hab : a ≤ b) (hbc : b ≤ c)
    (g : NeutralLogHilbertCarrier a) :
    ivIncl b c hbc (ivIncl a b hab g) = ivIncl a c (hab.trans hbc) g :=
  Subtype.ext rfl

theorem reconIncl_comp (a b c : ℝ) (hab : a ≤ b) (hbc : b ≤ c)
    (g : V08I.logCarrier a) :
    reconIncl b c hbc (reconIncl a b hab g) = reconIncl a c (hab.trans hbc) g :=
  Subtype.ext rfl

end WeilDefect.MasterReturn.LogCarrierNatural

#print axioms WeilDefect.MasterReturn.LogCarrierNatural.iv_mem_of_le
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.recon_mem_of_le
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.ivIncl
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.reconIncl
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.ivIncl_norm
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.reconIncl_norm
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.ivIncl_physical
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.reconIncl_physical
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.carrier_natural
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.ivIncl_comp
#print axioms WeilDefect.MasterReturn.LogCarrierNatural.reconIncl_comp
