import Mathlib

namespace CRALP.AnalyticGuard
set_option maxRecDepth 20000
set_option maxHeartbeats 0

/-- Transfer an exact natural cross-product inequality to a real remainder budget. -/
theorem ratio_lt_unit (N D S : ℕ) (hd : 0 < D) (hs : 0 < S)
    (h : N * S < D) : (N : ℝ) / (D : ℝ) < 1 / (S : ℝ) := by
  have hd' : (0 : ℝ) < D := by exact_mod_cast hd
  have hs' : (0 : ℝ) < S := by exact_mod_cast hs
  have h' : (N : ℝ) * (S : ℝ) < (D : ℝ) := by exact_mod_cast h
  apply (div_lt_div_iff₀ hd' hs').2
  simpa only [one_mul] using h'

theorem log_400_remainder : (2 * 7^801 * 400 : ℕ) * 10^340 < (801 * 351 * 20^801 : ℕ) := by decide

theorem log_400_remainder_real :
    ((2 * 7^801 * 400 : ℕ) : ℝ) / ((801 * 351 * 20^801 : ℕ) : ℝ) < 1 / ((10^340 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) log_400_remainder

theorem atan_400_remainder : (1 : ℕ) * 10^340 < (801 * 3^801 : ℕ) := by decide

theorem atan_400_remainder_real :
    ((1 : ℕ) : ℝ) / ((801 * 3^801 : ℕ) : ℝ) < 1 / ((10^340 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) atan_400_remainder

theorem sincos_400_remainder : (4^800 : ℕ) * 10^340 < (Nat.factorial 800 : ℕ) := by decide

theorem sincos_400_remainder_real :
    ((4^800 : ℕ) : ℝ) / ((Nat.factorial 800 : ℕ) : ℝ) < 1 / ((10^340 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) sincos_400_remainder

theorem exp_400_remainder : (3^400 * 401 : ℕ) * 10^340 < (398 * Nat.factorial 400 : ℕ) := by decide

theorem exp_400_remainder_real :
    ((3^400 * 401 : ℕ) : ℝ) / ((398 * Nat.factorial 400 : ℕ) : ℝ) < 1 / ((10^340 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) exp_400_remainder

theorem digamma_EM_512 : (4 * Nat.factorial 180 : ℕ) * 10^290 < (6^180 * 180 * 512^180 : ℕ) := by decide

theorem digamma_EM_512_real :
    ((4 * Nat.factorial 180 : ℕ) : ℝ) / ((6^180 * 180 * 512^180 : ℕ) : ℝ) < 1 / ((10^290 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) digamma_EM_512

theorem trigamma_EM_512 : (4 * Nat.factorial 180 : ℕ) * 10^290 < (6^180 * 512^181 : ℕ) := by decide

theorem trigamma_EM_512_real :
    ((4 * Nat.factorial 180 : ℕ) : ℝ) / ((6^180 * 512^181 : ℕ) : ℝ) < 1 / ((10^290 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) trigamma_EM_512

theorem digamma_EM_768 : (4 * Nat.factorial 180 : ℕ) * 10^290 < (6^180 * 180 * 768^180 : ℕ) := by decide

theorem digamma_EM_768_real :
    ((4 * Nat.factorial 180 : ℕ) : ℝ) / ((6^180 * 180 * 768^180 : ℕ) : ℝ) < 1 / ((10^290 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) digamma_EM_768

theorem trigamma_EM_768 : (4 * Nat.factorial 180 : ℕ) * 10^290 < (6^180 * 768^181 : ℕ) := by decide

theorem trigamma_EM_768_real :
    ((4 * Nat.factorial 180 : ℕ) : ℝ) / ((6^180 * 768^181 : ℕ) : ℝ) < 1 / ((10^290 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) trigamma_EM_768

theorem reflected_resolvent_K256_remainder : (100 * 2^1025 * 625 : ℕ) * 10^340 < (5^1025 * 609 : ℕ) := by decide

theorem reflected_resolvent_K256_remainder_real :
    ((100 * 2^1025 * 625 : ℕ) : ℝ) / ((5^1025 * 609 : ℕ) : ℝ) < 1 / ((10^340 : ℕ) : ℝ) := by
  exact ratio_lt_unit _ _ _ (by decide) (by decide) reflected_resolvent_K256_remainder

def expNum6 : ℕ := (11^0 * 10^6 * Nat.factorial 6 / Nat.factorial 0) + (11^1 * 10^5 * Nat.factorial 6 / Nat.factorial 1) + (11^2 * 10^4 * Nat.factorial 6 / Nat.factorial 2) + (11^3 * 10^3 * Nat.factorial 6 / Nat.factorial 3) + (11^4 * 10^2 * Nat.factorial 6 / Nat.factorial 4) + (11^5 * 10^1 * Nat.factorial 6 / Nat.factorial 5) + (11^6 * 10^0 * Nat.factorial 6 / Nat.factorial 6)
def expDen6 : ℕ := 10^6 * Nat.factorial 6
def expNum12 : ℕ := (11^0 * 20^12 * Nat.factorial 12 / Nat.factorial 0) + (11^1 * 20^11 * Nat.factorial 12 / Nat.factorial 1) + (11^2 * 20^10 * Nat.factorial 12 / Nat.factorial 2) + (11^3 * 20^9 * Nat.factorial 12 / Nat.factorial 3) + (11^4 * 20^8 * Nat.factorial 12 / Nat.factorial 4) + (11^5 * 20^7 * Nat.factorial 12 / Nat.factorial 5) + (11^6 * 20^6 * Nat.factorial 12 / Nat.factorial 6) + (11^7 * 20^5 * Nat.factorial 12 / Nat.factorial 7) + (11^8 * 20^4 * Nat.factorial 12 / Nat.factorial 8) + (11^9 * 20^3 * Nat.factorial 12 / Nat.factorial 9) + (11^10 * 20^2 * Nat.factorial 12 / Nat.factorial 10) + (11^11 * 20^1 * Nat.factorial 12 / Nat.factorial 11) + (11^12 * 20^0 * Nat.factorial 12 / Nat.factorial 12)
def expDen12 : ℕ := 20^12 * Nat.factorial 12
def expTailNum : ℕ := 11^13 * 280
def expTailDen : ℕ := 20^13 * Nat.factorial 13 * 269

theorem exp_minus_a_lt_2_over_5 : 5 * expDen6 < 2 * expNum6 := by decide
theorem exp_half_a_gt_17_over_10 : 17 * expDen12 < 10 * expNum12 := by decide
theorem exp_half_a_lt_7_over_4 : 4 * (expNum12 * expTailDen + expTailNum * expDen12) < 7 * expDen12 * expTailDen := by decide
theorem anchor_cosh_bound : (65 * 5 : ℕ) < 56 * 6 := by decide
theorem anchor_sinh_bound : (33 * 5 : ℕ) < 56 * 3 := by decide
theorem reflected_factor_uniform_bound : (16 : ℕ) < 100 := by decide

end CRALP.AnalyticGuard
#print axioms CRALP.AnalyticGuard.ratio_lt_unit
#print axioms CRALP.AnalyticGuard.log_400_remainder
#print axioms CRALP.AnalyticGuard.log_400_remainder_real
#print axioms CRALP.AnalyticGuard.atan_400_remainder
#print axioms CRALP.AnalyticGuard.atan_400_remainder_real
#print axioms CRALP.AnalyticGuard.sincos_400_remainder
#print axioms CRALP.AnalyticGuard.sincos_400_remainder_real
#print axioms CRALP.AnalyticGuard.exp_400_remainder
#print axioms CRALP.AnalyticGuard.exp_400_remainder_real
#print axioms CRALP.AnalyticGuard.digamma_EM_512
#print axioms CRALP.AnalyticGuard.digamma_EM_512_real
#print axioms CRALP.AnalyticGuard.trigamma_EM_512
#print axioms CRALP.AnalyticGuard.trigamma_EM_512_real
#print axioms CRALP.AnalyticGuard.digamma_EM_768
#print axioms CRALP.AnalyticGuard.digamma_EM_768_real
#print axioms CRALP.AnalyticGuard.trigamma_EM_768
#print axioms CRALP.AnalyticGuard.trigamma_EM_768_real
#print axioms CRALP.AnalyticGuard.reflected_resolvent_K256_remainder
#print axioms CRALP.AnalyticGuard.reflected_resolvent_K256_remainder_real
#print axioms CRALP.AnalyticGuard.exp_minus_a_lt_2_over_5
#print axioms CRALP.AnalyticGuard.exp_half_a_gt_17_over_10
#print axioms CRALP.AnalyticGuard.exp_half_a_lt_7_over_4
#print axioms CRALP.AnalyticGuard.anchor_cosh_bound
#print axioms CRALP.AnalyticGuard.anchor_sinh_bound
#print axioms CRALP.AnalyticGuard.reflected_factor_uniform_bound
