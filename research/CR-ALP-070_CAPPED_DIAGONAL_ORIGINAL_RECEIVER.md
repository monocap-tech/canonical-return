# CR-ALP-070 — capped diagonal original-form receiver

DRAFT / NONPROMOTING. CR-ALP-069 independent Fraction replay completed successfully (run 45625224037, job 139243501197): five original weighted artifact SHA256s match, independent inverse upper 0.06358596184791465 < K_lower 0.084343420037917, margin 0.02075745819000236. Native stdout recovered and read back at bb6c0b06c70a8dd35d9b0aba0cf6d9c388057ac0. This was the frozen old residual test, not original-form receiving acceptance.

## 1. Receiver preserving original correction

Let U=X'z be the original frozen even correction, r_L=Rz-LU and K_L the original corrected energy. Write C=L+A. Use the full 068 band diagonal rounded down to 1e-9, capped everywhere by its 2049 band value Dcap. This is a valid lower comparison and still dominates gamma on the entire Q85. Its uniform Weil floor dmin and tail floor dtail are unchanged from 069. A=Q85 diag(D-gamma) Q85 is bounded positive with norm <= Acap=Dcap-gamma. Correct receiving energy and residual are K_C=K_L+<U,AU>, r_C=r_L-AU. These apply to the original one-dimensional frozen source direction, not the full source bank.

The archived 058 contracted correction bank contains 425 even sine coefficients t_j and a pole coefficient db, in input_2.txt. Define U0=Q85 sum t_j e_j using the exact physical projection. Let st=||b_high||², c_j=<b,e_j>, tau=sum_(j>85)t_j c_j/st. Then U0_j=t_j-tau c_j for j>85, with t_j=0 beyond 850. Check the archive db overlaps -tau; this is a consistency check, not a proof of original correction attachment. For j>2048 all A_j=Acap. Consequently AU0=Q85 z0, where z0 is supported on 86..2048 and

  z0_j=A_j t_j-(A_j-Acap)tau c_j.

Indeed subtracting the entire constant multiple of b_high from the unprojected A action is annihilated by Q85. Thus <U0,AU0>=sum_(86<=j<=2048) U0_j z0_j exactly. No high-sine truncation of the original correction is made.

## 2. Algebraic response transfer including original receiver

Let x be the recovered 067 weighted trial. Delta=(D_new-D_old)x. For its archived physical response residual y=rho-(D_old+P_rat)x and energy E_old, update y_new=y-Delta-z0, E_new=E_old-<x,Delta>-2<x,z0>. Q85 y_new is the response residual for rho-AU0. Compute full norm and pole pairing from existing logs by finite algebraic cross terms. Complete tail squared norm comes from Parseval after exact Q85 projection, not a finite cutoff assumption. The energy and full weighted inverse residual upper bound are computed using Fraction endpoints. Retain the 1e-13 shift displacement allowance and original proxy L2 error.

The native test explicitly assumes an attachment tolerance deltaX=1e-10 for ||U-U0||. This is a proposed sufficient tolerance, not an established premise. The source attachment error adds Acap*deltaX to the inverse residual proxy error, and subtracts at most Acap*deltaX*(2||U0||+deltaX) from the lower energy correction. Since A>=0, the lower energy correction is at least zero independently. The archive finite high sine norm bounds ||U0|| by projection contractivity. If the test passes, it establishes a conditional original-form scalar receiver provided this source attachment is discharged. It does not establish a full matrix Schur certificate, new index decrement, or RH.

## 3. Implementation and custody

src/audits/cr_alp_070_receiver.py committed/read back at 0db04a3b8d76c91c8b7b58b8e03081427faeafc7. It authenticates five original weighted artifacts, reconstructs exact original response trial from its rational bank, reads authenticated archive input_2.txt, checks pole bridge overlap, and evaluates the full conditional receiver with outward Fraction endpoints. No additional physical integration or candidate solve. Workflow .github/workflows/cr-alp-070.yml committed/read back at c4d049a21caba76bb14b21a8f2cc95661398b9be. It reassembles and checks original archive, prints exact preparation/assembly source and original correction bank for source-attachment review, and runs the candidate. Native execution not yet verified. The original correction tolerance must be derived from archived generation, not inferred from the pole overlap or residual proxy error alone. Preserve source bank/logs and output before further promotion.

Prior odd frozen scalar comparison, 1048-dimensional positivity and conditional nonpositive-index bound <=66 unchanged. Historical standalone archive reconciliation remains open. PR #3 stays draft.
