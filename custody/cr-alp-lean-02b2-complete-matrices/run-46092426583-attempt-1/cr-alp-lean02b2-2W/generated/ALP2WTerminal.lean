import ALP2WRows00
import ALP2WRows05
import ALP2WRows10
import ALP2WRows15
import ALP2WRows20
import ALP2WRows25
import ALP2WRows30
import ALP2WRows35
import ALP2WRows40
namespace CRALP.Cert2W
open scoped BigOperators
set_option maxRecDepth 20000
set_option maxHeartbeats 0
theorem allRows : ∀ i, rowOK M N D S d p i := by
  intro i
  fin_cases i
  · exact row_0
  · exact row_1
  · exact row_2
  · exact row_3
  · exact row_4
  · exact row_5
  · exact row_6
  · exact row_7
  · exact row_8
  · exact row_9
  · exact row_10
  · exact row_11
  · exact row_12
  · exact row_13
  · exact row_14
  · exact row_15
  · exact row_16
  · exact row_17
  · exact row_18
  · exact row_19
  · exact row_20
  · exact row_21
  · exact row_22
  · exact row_23
  · exact row_24
  · exact row_25
  · exact row_26
  · exact row_27
  · exact row_28
  · exact row_29
  · exact row_30
  · exact row_31
  · exact row_32
  · exact row_33
  · exact row_34
  · exact row_35
  · exact row_36
  · exact row_37
  · exact row_38
  · exact row_39
  · exact row_40
  · exact row_41
theorem matrixLower (x : Fin 42 → ℝ) :
    ((p : ℝ) / (2 * d)) * (∑ i, x i ^ 2) ≤
      quad (fun i j => (M i j : ℝ) / D) x := by
  exact intWitnessLower M N D S d p (by decide) (by decide) (by decide)
    (by decide) (by decide) allRows x
end CRALP.Cert2W
#print axioms CRALP.Cert2W.allRows
#print axioms CRALP.Cert2W.matrixLower
