import Mathlib

namespace CRALP.DirectedInterval

/-- Real containment; bounds are the actual scaled endpoint values. -/
def Contains (lo hi x : ℝ) : Prop := lo ≤ x ∧ x ≤ hi

theorem widen {a b x l u : ℝ} (h : Contains a b x)
    (hl : l ≤ a) (hu : b ≤ u) : Contains l u x :=
  ⟨hl.trans h.1, h.2.trans hu⟩

theorem add {a b c d x y : ℝ} (hx : Contains a b x)
    (hy : Contains c d y) : Contains (a+c) (b+d) (x+y) :=
  ⟨add_le_add hx.1 hy.1, add_le_add hx.2 hy.2⟩

theorem neg {a b x : ℝ} (hx : Contains a b x) :
    Contains (-b) (-a) (-x) :=
  ⟨neg_le_neg hx.2, neg_le_neg hx.1⟩

theorem sub {a b c d x y : ℝ} (hx : Contains a b x)
    (hy : Contains c d y) : Contains (a-d) (b-c) (x-y) := by
  simpa only [sub_eq_add_neg] using add hx (neg hy)

/-- Both endpoints of a linear function bound every interior value. -/
theorem scaleCorners {c d y t l u : ℝ} (hy : Contains c d y)
    (hc : Contains l u (t*c)) (hd : Contains l u (t*d)) :
    Contains l u (t*y) := by
  by_cases ht : 0 ≤ t
  · exact ⟨hc.1.trans (mul_le_mul_of_nonneg_left hy.1 ht),
      (mul_le_mul_of_nonneg_left hy.2 ht).trans hd.2⟩
  · have ht' : t ≤ 0 := le_of_not_ge ht
    exact ⟨hd.1.trans (mul_le_mul_of_nonpos_left hy.2 ht'),
      (mul_le_mul_of_nonpos_left hy.1 ht').trans hc.2⟩

/-- Four corners with outward rounding suffice for multiplication. -/
theorem mulCorners {a b c d x y l u : ℝ}
    (hx : Contains a b x) (hy : Contains c d y)
    (hac : Contains l u (a*c)) (had : Contains l u (a*d))
    (hbc : Contains l u (b*c)) (hbd : Contains l u (b*d)) :
    Contains l u (x*y) := by
  have ha := scaleCorners hy hac had
  have hb := scaleCorners hy hbc hbd
  by_cases h : 0 ≤ y
  · exact ⟨ha.1.trans (mul_le_mul_of_nonneg_right hx.1 h),
      (mul_le_mul_of_nonneg_right hx.2 h).trans hb.2⟩
  · have h' : y ≤ 0 := le_of_not_ge h
    exact ⟨hb.1.trans (mul_le_mul_of_nonpos_right hx.2 h'),
      (mul_le_mul_of_nonpos_right hx.1 h').trans ha.2⟩

theorem reciprocalPositive {c d y : ℝ} (hc : 0 < c)
    (hy : Contains c d y) : Contains (1/d) (1/c) (1/y) := by
  have hy' : 0 < y := lt_of_lt_of_le hc hy.1
  have hd : 0 < d := lt_of_lt_of_le hy' hy.2
  constructor
  · apply (div_le_div_iff₀ hd hy').2
    simpa only [one_mul] using hy.2
  · apply (div_le_div_iff₀ hy' hc).2
    simpa only [one_mul] using hy.1

theorem divCornersPositive {a b c d x y l u : ℝ}
    (hx : Contains a b x) (hy : Contains c d y) (hc : 0 < c)
    (hac : Contains l u (a/c)) (had : Contains l u (a/d))
    (hbc : Contains l u (b/c)) (hbd : Contains l u (b/d)) :
    Contains l u (x/y) := by
  have h := mulCorners hx (reciprocalPositive hc hy)
    (show Contains l u (a*(1/d)) by simpa only [mul_one_div] using had)
    (show Contains l u (a*(1/c)) by simpa only [mul_one_div] using hac)
    (show Contains l u (b*(1/d)) by simpa only [mul_one_div] using hbd)
    (show Contains l u (b*(1/c)) by simpa only [mul_one_div] using hbc)
  simpa only [mul_one_div] using h

theorem divCornersNegative {a b c d x y l u : ℝ}
    (hx : Contains a b x) (hy : Contains c d y) (hd : d < 0)
    (hac : Contains l u (a/c)) (had : Contains l u (a/d))
    (hbc : Contains l u (b/c)) (hbd : Contains l u (b/d)) :
    Contains l u (x/y) := by
  have h := divCornersPositive (neg hx) (neg hy) (neg_pos.mpr hd)
    (show Contains l u ((-b)/(-d)) by simpa only [neg_div_neg_eq] using hbd)
    (show Contains l u ((-b)/(-c)) by simpa only [neg_div_neg_eq] using hbc)
    (show Contains l u ((-a)/(-d)) by simpa only [neg_div_neg_eq] using had)
    (show Contains l u ((-a)/(-c)) by simpa only [neg_div_neg_eq] using hac)
  simpa only [neg_div_neg_eq] using h

theorem squarePositive {a b x : ℝ} (ha : 0 ≤ a)
    (hx : Contains a b x) : Contains (a^2) (b^2) (x^2) := by
  have hx0 := ha.trans hx.1
  have hb0 := hx0.trans hx.2
  constructor
  · nlinarith [mul_nonneg (sub_nonneg.mpr hx.1) (show 0 ≤ x+a by linarith)]
  · nlinarith [mul_nonneg (sub_nonneg.mpr hx.2) (show 0 ≤ b+x by linarith)]

theorem squareNegative {a b x : ℝ} (hb : b ≤ 0)
    (hx : Contains a b x) : Contains (b^2) (a^2) (x^2) := by
  have h := squarePositive (neg_nonneg.mpr hb) (neg hx)
  simpa only [neg_sq] using h

theorem squareCrossZero {a b x : ℝ} (ha : a ≤ 0) (hb : 0 ≤ b)
    (hx : Contains a b x) : Contains 0 (max (a^2) (b^2)) (x^2) := by
  refine ⟨sq_nonneg x, ?_⟩
  by_cases hp : 0 ≤ x
  · exact (squarePositive (le_refl 0) (show Contains 0 b x from ⟨hp,hx.2⟩)).2.trans (le_max_right _ _)
  · exact (squareNegative (le_refl 0) (show Contains a 0 x from ⟨hx.1,le_of_not_ge hp⟩)).2.trans (le_max_left _ _)

theorem sqrtRounded {a b x l u : ℝ} (hx : Contains a b x)
    (hl : 0 ≤ l) (hu : 0 ≤ u) (hlo : l^2 ≤ a) (hhi : b ≤ u^2) :
    Contains l u (Real.sqrt x) := by
  constructor
  · calc
      l = Real.sqrt (l^2) := by rw [Real.sqrt_sq_eq_abs, abs_of_nonneg hl]
      _ ≤ Real.sqrt x := Real.sqrt_le_sqrt (hlo.trans hx.1)
  · calc
      Real.sqrt x ≤ Real.sqrt (u^2) := Real.sqrt_le_sqrt (hx.2.trans hhi)
      _ = u := by rw [Real.sqrt_sq_eq_abs, abs_of_nonneg hu]

/-- An exact signed cross-product certificate gives directed endpoint inclusion.
    This is a receiver for integer floor/ceiling witnesses, not a claim
    that any particular external rounding trace has been verified. -/
theorem scaledRatio (l u n d s : ℤ) (hd : 0 < d) (hs : 0 < s)
    (hl : l*d ≤ n*s) (hu : n*s ≤ u*d) :
    Contains ((l : ℝ)/s) ((u : ℝ)/s) ((n : ℝ)/d) := by
  have hd' : (0 : ℝ) < d := by exact_mod_cast hd
  have hs' : (0 : ℝ) < s := by exact_mod_cast hs
  constructor
  · apply (div_le_div_iff₀ hs' hd').2
    exact_mod_cast hl
  · apply (div_le_div_iff₀ hd' hs').2
    exact_mod_cast hu

end CRALP.DirectedInterval
#print axioms CRALP.DirectedInterval.widen
#print axioms CRALP.DirectedInterval.add
#print axioms CRALP.DirectedInterval.neg
#print axioms CRALP.DirectedInterval.sub
#print axioms CRALP.DirectedInterval.scaleCorners
#print axioms CRALP.DirectedInterval.mulCorners
#print axioms CRALP.DirectedInterval.reciprocalPositive
#print axioms CRALP.DirectedInterval.divCornersPositive
#print axioms CRALP.DirectedInterval.divCornersNegative
#print axioms CRALP.DirectedInterval.squarePositive
#print axioms CRALP.DirectedInterval.squareNegative
#print axioms CRALP.DirectedInterval.squareCrossZero
#print axioms CRALP.DirectedInterval.sqrtRounded
#print axioms CRALP.DirectedInterval.scaledRatio
