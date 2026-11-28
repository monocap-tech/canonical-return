import GridTrace

namespace CRALP.SignedGridSeries
open Finset
open scoped BigOperators
open CRALP.GridTrace CRALP.DirectedEvaluation CRALP.DirectedInterval CRALP.RationalTrace

theorem sumSound (t : ℕ → Grid) (f : ℕ → ℝ)
    (h : ∀ i ∈ range 400, Encl (t i) (f i)) :
    Contains (val (∑ i ∈ range 400, (t i).lo))
      (val (∑ i ∈ range 400, (t i).hi)) (∑ i ∈ range 400, f i) := by
  have hl : (∑ i ∈ range 400, val (t i).lo) ≤ ∑ i ∈ range 400, f i :=
    Finset.sum_le_sum (fun i hi => (h i hi).1)
  have hu : (∑ i ∈ range 400, f i) ≤ ∑ i ∈ range 400, val (t i).hi :=
    Finset.sum_le_sum (fun i hi => (h i hi).2)
  unfold val at *
  push_cast
  rw [Finset.sum_div,Finset.sum_div]
  exact ⟨hl,hu⟩

theorem attachGrid {t : ℕ → Grid} {f : ℕ → ℝ} {y : ℝ} {out : Grid}
    (h : ∀ i ∈ range 400, Encl (t i) (f i))
    (he : |y-∑ i ∈ range 400, f i| < ElementarySeries.unit)
    (hl : out.lo=(∑ i ∈ range 400, (t i).lo)-1)
    (hu : out.hi=(∑ i ∈ range 400, (t i).hi)+1) :
    Encl out y := by
  have ha := attach (sumSound t f h) he
  have hL : val out.lo = val (∑ i ∈ range 400, (t i).lo)-ElementarySeries.unit := by
    rw [hl,unitEq]
    unfold val
    push_cast <;> ring
  have hU : val out.hi = val (∑ i ∈ range 400, (t i).hi)+ElementarySeries.unit := by
    rw [hu,unitEq]
    unfold val
    push_cast <;> ring
  unfold Encl
  rw [hL,hU]
  exact ha

theorem atanPowerTrace {x : ℝ} (v : Grid) (t term : ℕ → Grid)
    (ht : Encl (t 0) x) (hv : Encl v (-(x^2)))
    (hstep : ∀ i : Fin 400, logStep v t term i) :
    ∀ i ≤ 400, Encl (t i) ((-1:ℝ)^i*x^(2*i+1)) := by
  intro i
  induction i with
  | zero =>
    intro _
    simpa using ht
  | succ i ih =>
    intro hi
    have hj : i < 400 := by omega
    have h := mulSound (ih (by omega)) hv (hstep ⟨i,hj⟩).2
    rw [pow_succ (-1:ℝ) i,show 2*(i+1)+1=2*i+1+2 by omega,pow_add]
    convert h using 1 <;> ring

theorem atanTrace {x : ℝ} (v : Grid) (t term : ℕ → Grid) (out : Grid)
    (hr : |x| ≤ 1/3) (ht : Encl (t 0) x) (hv : Encl v (-(x^2)))
    (hstep : ∀ i : Fin 400, logStep v t term i)
    (hl : out.lo=(∑ i ∈ range 400, (term i).lo)-1)
    (hu : out.hi=(∑ i ∈ range 400, (term i).hi)+1) :
    Encl out (Real.arctan x) := by
  have hp := atanPowerTrace v t term ht hv hstep
  have hterm (i : ℕ) (hi : i ∈ range 400) :
      Encl (term i) ((-1:ℝ)^i*x^(2*i+1)/(2*i+1 : ℕ)) := by
    have hj := mem_range.mp hi
    have h := divSound (hp i (by omega)) (pointSound ((2*i+1 : ℕ):ℤ)) (hstep ⟨i,hj⟩).1
    simpa only [Int.cast_natCast] using h
  exact attachGrid hterm (ElementarySeries.atanReduced x hr) hl hu

noncomputable def trigTerm (e i : ℕ) (x : ℝ) : ℝ :=
  (-1:ℝ)^i*x^(2*i+e)/(Nat.factorial (2*i+e):ℝ)

def trigDen (e i : ℕ) : ℕ := (2*i+e+1)*(2*i+e+2)

theorem trigTermSucc (e i : ℕ) (x : ℝ) :
    trigTerm e (i+1) x = trigTerm e i x * (-(x^2)) / (trigDen e i:ℝ) := by
  unfold trigTerm trigDen
  rw [show 2*(i+1)+e=2*i+e+1+1 by omega,
      Nat.factorial_succ,Nat.factorial_succ]
  push_cast
  rw [show 2*i+e+1+1=2*i+e+2 by omega,pow_add,pow_succ (-1:ℝ) i]
  have hf : (Nat.factorial (2*i+e):ℝ) ≠ 0 := by positivity
  have h1 : (2*(i:ℝ)+(e:ℝ)+1) ≠ 0 := by positivity
  have h2 : (2*(i:ℝ)+(e:ℝ)+1+1) ≠ 0 := by positivity
  have h3 : (2*(i:ℝ)+(e:ℝ)+2) ≠ 0 := by positivity
  field_simp [hf,h1,h2,h3] <;> ring

def trigStep (e : ℕ) (v : Grid) (t mid : ℕ → Grid) (i : Fin 399) : Prop :=
  mulOK (t i.val) v (mid i.val) ∧
  divOK (mid i.val) (point (trigDen e i.val:ℤ)) (t (i.val+1))

instance (e : ℕ) (v : Grid) (t mid : ℕ → Grid) (i : Fin 399) :
    Decidable (trigStep e v t mid i) := by
  unfold trigStep
  infer_instance

theorem trigPowerTrace {x : ℝ} (e : ℕ) (v : Grid) (t mid : ℕ → Grid)
    (ht : Encl (t 0) (trigTerm e 0 x)) (hv : Encl v (-(x^2)))
    (hstep : ∀ i : Fin 399, trigStep e v t mid i) :
    ∀ i < 400, Encl (t i) (trigTerm e i x) := by
  intro i
  induction i with
  | zero => intro _; exact ht
  | succ i ih =>
    intro hi
    have hj : i < 399 := by omega
    have hs := hstep ⟨i,hj⟩
    have hm := mulSound (ih (by omega)) hv hs.1
    have h := divSound hm (pointSound (trigDen e i:ℤ)) hs.2
    rw [trigTermSucc]
    simpa only [Int.cast_natCast] using h

theorem sinTrace {x : ℝ} (v : Grid) (t mid : ℕ → Grid) (out : Grid)
    (hr : |x| ≤ 4) (ht : Encl (t 0) x) (hv : Encl v (-(x^2)))
    (hstep : ∀ i : Fin 399, trigStep 1 v t mid i)
    (hl : out.lo=(∑ i ∈ range 400, (t i).lo)-1)
    (hu : out.hi=(∑ i ∈ range 400, (t i).hi)+1) :
    Encl out (Real.sin x) := by
  have ht' : Encl (t 0) (trigTerm 1 0 x) := by simpa [trigTerm] using ht
  have hp := trigPowerTrace 1 v t mid ht' hv hstep
  apply attachGrid (fun i hi => hp i (mem_range.mp hi)) _ hl hu
  have h := TrigSeries.sinReduced x hr
  simpa only [TrigSeries.trigUnit,ElementarySeries.unit,TrigSeries.sinPoly,trigTerm] using h

theorem cosTrace {x : ℝ} (v : Grid) (t mid : ℕ → Grid) (out : Grid)
    (hr : |x| ≤ 4) (ht : Encl (t 0) 1) (hv : Encl v (-(x^2)))
    (hstep : ∀ i : Fin 399, trigStep 0 v t mid i)
    (hl : out.lo=(∑ i ∈ range 400, (t i).lo)-1)
    (hu : out.hi=(∑ i ∈ range 400, (t i).hi)+1) :
    Encl out (Real.cos x) := by
  have ht' : Encl (t 0) (trigTerm 0 0 x) := by simpa [trigTerm] using ht
  have hp := trigPowerTrace 0 v t mid ht' hv hstep
  apply attachGrid (fun i hi => hp i (mem_range.mp hi)) _ hl hu
  have h := TrigSeries.cosReduced x hr
  simpa only [TrigSeries.trigUnit,ElementarySeries.unit,TrigSeries.cosPoly,trigTerm,Nat.add_zero] using h

end CRALP.SignedGridSeries
#print axioms CRALP.SignedGridSeries.sumSound
#print axioms CRALP.SignedGridSeries.attachGrid
#print axioms CRALP.SignedGridSeries.atanPowerTrace
#print axioms CRALP.SignedGridSeries.atanTrace
#print axioms CRALP.SignedGridSeries.trigTermSucc
#print axioms CRALP.SignedGridSeries.trigPowerTrace
#print axioms CRALP.SignedGridSeries.sinTrace
#print axioms CRALP.SignedGridSeries.cosTrace
