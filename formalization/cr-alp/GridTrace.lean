import RationalTrace

namespace CRALP.GridTrace
open Finset
open scoped BigOperators
open CRALP.DirectedInterval CRALP.DirectedEvaluation CRALP.RationalTrace

structure Grid where
  lo : ℤ
  hi : ℤ
  deriving DecidableEq

def val (z : ℤ) : ℝ := (z:ℝ)/logScale
def Encl (a : Grid) (x : ℝ) : Prop := Contains (val a.lo) (val a.hi) x

theorem scalePos : 0 < logScale := by unfold logScale; decide
theorem scalePosReal : (0:ℝ) < logScale := by exact_mod_cast scalePos

def point (z : ℤ) : Grid := ⟨z*logScale,z*logScale⟩
def plus (a b : Grid) : Grid := ⟨a.lo+b.lo,a.hi+b.hi⟩
def minus (a b : Grid) : Grid := ⟨a.lo-b.hi,a.hi-b.lo⟩
def negative (a : Grid) : Grid := ⟨-a.hi,-a.lo⟩

theorem pointSound (z : ℤ) : Encl (point z) (z:ℝ) := by
  have he : val (z*logScale) = (z:ℝ) := by
    unfold val
    push_cast
    field_simp [scalePosReal.ne']
  exact ⟨by simpa only [point,he],by simpa only [point,he]⟩

theorem plusSound {a b : Grid} {x y : ℝ} (ha : Encl a x) (hb : Encl b y) :
    Encl (plus a b) (x+y) := by
  unfold Encl plus val at *
  push_cast
  rw [add_div,add_div]
  exact DirectedInterval.add ha hb

theorem minusSound {a b : Grid} {x y : ℝ} (ha : Encl a x) (hb : Encl b y) :
    Encl (minus a b) (x-y) := by
  unfold Encl minus val at *
  push_cast
  rw [sub_div,sub_div]
  exact DirectedInterval.sub ha hb

theorem negativeSound {a : Grid} {x : ℝ} (ha : Encl a x) :
    Encl (negative a) (-x) := by
  unfold Encl negative val at *
  push_cast
  rw [neg_div,neg_div]
  exact DirectedInterval.neg ha

def mulCornerOK (a b : ℤ) (c : Grid) : Prop :=
  c.lo*logScale ≤ a*b ∧ a*b ≤ c.hi*logScale

instance (a b : ℤ) (c : Grid) : Decidable (mulCornerOK a b c) := by
  unfold mulCornerOK
  infer_instance

theorem mulCornerSound {a b : ℤ} {c : Grid} (h : mulCornerOK a b c) :
    Encl c (val a*val b) := by
  have hl := mul_le_mul_of_nonneg_right h.1 scalePos.le
  have hu := mul_le_mul_of_nonneg_right h.2 scalePos.le
  have hc := scaledRatio c.lo c.hi (a*b) (logScale*logScale) logScale
    (mul_pos scalePos scalePos) scalePos
    (by simpa only [mul_assoc] using hl) (by simpa only [mul_assoc] using hu)
  unfold Encl val
  rw [div_mul_div_comm]
  simpa only [Int.cast_mul] using hc

def mulOK (a b c : Grid) : Prop :=
  mulCornerOK a.lo b.lo c ∧ mulCornerOK a.lo b.hi c ∧
  mulCornerOK a.hi b.lo c ∧ mulCornerOK a.hi b.hi c

instance (a b c : Grid) : Decidable (mulOK a b c) := by
  unfold mulOK
  infer_instance

theorem mulSound {a b c : Grid} {x y : ℝ}
    (ha : Encl a x) (hb : Encl b y) (h : mulOK a b c) : Encl c (x*y) :=
  mulCorners ha hb (mulCornerSound h.1) (mulCornerSound h.2.1)
    (mulCornerSound h.2.2.1) (mulCornerSound h.2.2.2)

def divCornerOK (a b : ℤ) (c : Grid) : Prop :=
  0 < b ∧ c.lo*b ≤ a*logScale ∧ a*logScale ≤ c.hi*b

instance (a b : ℤ) (c : Grid) : Decidable (divCornerOK a b c) := by
  unfold divCornerOK
  infer_instance

theorem divCornerSound {a b : ℤ} {c : Grid} (h : divCornerOK a b c) :
    Encl c (val a/val b) := by
  have hb : (0:ℝ) < b := by exact_mod_cast h.1
  have he : val a/val b = (a:ℝ)/(b:ℝ) := by
    unfold val
    field_simp [scalePosReal.ne',hb.ne'] <;> ring
  rw [he]
  exact scaledRatio c.lo c.hi a b logScale h.1 scalePos h.2.1 h.2.2

def divOK (a b c : Grid) : Prop :=
  divCornerOK a.lo b.lo c ∧ divCornerOK a.lo b.hi c ∧
  divCornerOK a.hi b.lo c ∧ divCornerOK a.hi b.hi c

instance (a b c : Grid) : Decidable (divOK a b c) := by
  unfold divOK
  infer_instance

theorem divSound {a b c : Grid} {x y : ℝ}
    (ha : Encl a x) (hb : Encl b y) (h : divOK a b c) : Encl c (x/y) := by
  have hb0 : (0:ℝ) < b.lo := by exact_mod_cast h.1.1
  exact divCornersPositive ha hb (div_pos hb0 scalePosReal)
    (divCornerSound h.1) (divCornerSound h.2.1)
    (divCornerSound h.2.2.1) (divCornerSound h.2.2.2)

def absGuard (a : Grid) (p q : ℤ) : Prop :=
  0 < q ∧ -p*logScale ≤ a.lo*q ∧ a.hi*q ≤ p*logScale

instance (a : Grid) (p q : ℤ) : Decidable (absGuard a p q) := by
  unfold absGuard
  infer_instance

theorem absSound {a : Grid} {x : ℝ} {p q : ℤ}
    (ha : Encl a x) (h : absGuard a p q) : |x| ≤ (p:ℝ)/q := by
  have hq : (0:ℝ) < q := by exact_mod_cast h.1
  have hl : -(p:ℝ)/(q:ℝ) ≤ val a.lo := by
    unfold val
    apply (div_le_div_iff₀ hq scalePosReal).2
    exact_mod_cast h.2.1
  have hu : val a.hi ≤ (p:ℝ)/(q:ℝ) := by
    unfold val
    apply (div_le_div_iff₀ scalePosReal hq).2
    exact_mod_cast h.2.2
  exact abs_le.2 ⟨by simpa only [neg_div] using hl.trans ha.1,ha.2.trans hu⟩

def logStep (v : Grid) (t term : ℕ → Grid) (i : Fin 400) : Prop :=
  divOK (t i.val) (point ((2*i.val+1 : ℕ):ℤ)) (term i.val) ∧
  mulOK (t i.val) v (t (i.val+1))

instance (v : Grid) (t term : ℕ → Grid) (i : Fin 400) : Decidable (logStep v t term i) := by
  unfold logStep
  infer_instance

theorem logPowerTrace {u : ℝ} (v : Grid) (t term : ℕ → Grid)
    (ht : Encl (t 0) u) (hv : Encl v (u^2))
    (hstep : ∀ i : Fin 400, logStep v t term i) :
    ∀ i ≤ 400, Encl (t i) (u^(2*i+1)) := by
  intro i
  induction i with
  | zero =>
    intro _
    simpa using ht
  | succ i ih =>
    intro hi
    have hj : i < 400 := by omega
    have h := mulSound (ih (by omega)) hv (hstep ⟨i,hj⟩).2
    rw [show 2*(i+1)+1=2*i+1+2 by omega,pow_add]
    exact h

theorem logPolynomialTrace {u : ℝ} (v : Grid) (t term : ℕ → Grid)
    (ht : Encl (t 0) u) (hv : Encl v (u^2))
    (hstep : ∀ i : Fin 400, logStep v t term i) :
    Contains (2*val (∑ i ∈ range 400, (term i).lo))
      (2*val (∑ i ∈ range 400, (term i).hi)) (ElementarySeries.logPoly u) := by
  have hp := logPowerTrace v t term ht hv hstep
  have hterm (i : ℕ) (hi : i ∈ range 400) :
      Encl (term i) (u^(2*i+1)/(2*i+1 : ℕ)) := by
    have hj := mem_range.mp hi
    have h := divSound (hp i (by omega)) (pointSound ((2*i+1 : ℕ):ℤ)) (hstep ⟨i,hj⟩).1
    simpa only [Int.cast_natCast] using h
  have hl : (∑ i ∈ range 400, val (term i).lo) ≤
      ∑ i ∈ range 400, u^(2*i+1)/(2*i+1 : ℕ) :=
    Finset.sum_le_sum (fun i hi => (hterm i hi).1)
  have hu : (∑ i ∈ range 400, u^(2*i+1)/(2*i+1 : ℕ)) ≤
      ∑ i ∈ range 400, val (term i).hi :=
    Finset.sum_le_sum (fun i hi => (hterm i hi).2)
  unfold ElementarySeries.logPoly val at *
  push_cast
  rw [Finset.sum_div,Finset.sum_div]
  exact ⟨mul_le_mul_of_nonneg_left hl (by norm_num),mul_le_mul_of_nonneg_left hu (by norm_num)⟩

theorem logTrace {x : ℝ} (v : Grid) (t term : ℕ → Grid) (out : Grid)
    (hx : 0 < x) (hr : |(x-1)/(x+1)| ≤ 7/20)
    (ht : Encl (t 0) ((x-1)/(x+1))) (hv : Encl v (((x-1)/(x+1))^2))
    (hstep : ∀ i : Fin 400, logStep v t term i)
    (hl : out.lo=2*(∑ i ∈ range 400, (term i).lo)-1)
    (hu : out.hi=2*(∑ i ∈ range 400, (term i).hi)+1) :
    Encl out (Real.log x) := by
  have h := logAttach hx hr (logPolynomialTrace v t term ht hv hstep)
  unfold Encl val
  rw [hl,hu]
  push_cast
  rw [sub_div,add_div,← unitEq]
  convert h using 1 <;> unfold val <;> push_cast <;> ring

end CRALP.GridTrace
#print axioms CRALP.GridTrace.scalePos
#print axioms CRALP.GridTrace.scalePosReal
#print axioms CRALP.GridTrace.pointSound
#print axioms CRALP.GridTrace.plusSound
#print axioms CRALP.GridTrace.minusSound
#print axioms CRALP.GridTrace.negativeSound
#print axioms CRALP.GridTrace.mulCornerSound
#print axioms CRALP.GridTrace.mulSound
#print axioms CRALP.GridTrace.divCornerSound
#print axioms CRALP.GridTrace.divSound
#print axioms CRALP.GridTrace.absSound
#print axioms CRALP.GridTrace.logPowerTrace
#print axioms CRALP.GridTrace.logPolynomialTrace
#print axioms CRALP.GridTrace.logTrace
