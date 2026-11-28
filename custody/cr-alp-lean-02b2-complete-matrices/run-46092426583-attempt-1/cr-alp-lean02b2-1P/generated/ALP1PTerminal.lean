import ALP1PRows00
import ALP1PRows05
import ALP1PRows10
import ALP1PRows15
import ALP1PRows20
import ALP1PRows25
import ALP1PRows30
import ALP1PRows35
import ALP1PRows40
namespace CRALP.Cert1P
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
  · exact row_42
theorem matrixLower (x : Fin 43 → ℝ) :
    ((p : ℝ) / (2 * d)) * (∑ i, x i ^ 2) ≤
      quad (fun i j => (M i j : ℝ) / D) x := by
  exact intWitnessLower M N D S d p (by decide) (by decide) (by decide)
    (by decide) (by decide) allRows x
end CRALP.Cert1P
#print axioms CRALP.Cert1P.allRows
#print axioms CRALP.Cert1P.matrixLower
