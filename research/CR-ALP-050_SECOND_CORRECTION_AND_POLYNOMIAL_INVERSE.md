# CR-ALP-050 — Second structured correction, 408 finite directions, and Chebyshev inverse residual

DRAFT / NONPROMOTING. Exact local arithmetic conditional on inherited continuum identities. No complete infinite Schur, Lean, or RH claim.

## Inputs and execution

Reuse CR-ALP-046 E85 bank SHA256 d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22 and CR-ALP-049 saved rational V and Y, whose input SHA256s are pinned in audit.py. CR-ALP-045 layered Gamma lower and CR-ALP-042 signed-prime norm imply L=gamma I+Q85 P Q85 with gamma=72/25, a0=51/200 <= L <= b0=1101/200. The structured L is the same as CR-ALP-048/049. No automatic transfer from NF662.

Use matching-parity ordinary high sine rows 86..600, and the exact pole-null trial through 410. Embed the previous 205 trial X0 by zero extension in the enlarged ordinary sine band. Solve the finite constrained saddle system [Ltt c; c^T 0] for the correction of U-Ltt X0 using double precision only to choose coordinates. Round free rows to multiples of 10^-16, impose the last row by the genuine pole coefficient ratios. X=X0+DX is exactly pole-null (the last row is an interval enclosing the exact analytic coefficient). Preserve the same physical low V as CR-ALP-049; exact interval energy congruence remains within 10^-12 of I.

The inherited signed coefficient generator uses 100 decimal outward integers, 24 Bernoulli terms and shift 1024 (increased from 512 to support index 600); analytic digamma/trigamma remainder is bounded as in the code. All seven signed prime translations are retained. For exact interval matrix multiplication, use integer midpoint/radius products M_A M_B and |M_A|R_B+R_A|M_B|+R_A R_B, then outward division by the fixed-point scale. This is an entrywise enclosure, not a floating point acceptance test. The executed generator and result banks are in this directory; the extra predecessor hash assertions were inserted after the run and their input hashes checked separately.

## 408-dimensional finite structured comparison

Set K=V^T H V-U^T X-X^T U+X^T Ltt X, where U is the ordinary signed coefficient matrix through the trial band, and X is supported there. The complete in-band pole projection gives G_trial=Rtrial^T Rtrial-(Rtrial^T c)(c^T Rtrial)/(c^T c), where Rtrial=U-Ltt X. Both parities have ||G_trial||<10^-26 by absolute row sums. Exact shifted midpoint/radius Bareiss tests give K>(1/4)I in odd-index parity (43 pivots), K>(1/8)I in even-index parity (42 pivots). Diagnostic midpoint minimum eigenvalues of the new K are 0.2936857504 and 0.2013114650, not certified eigenvalue enclosures.

The restricted structured operator L_F>=51/200 I. Its exact Schur complement equals K-Rtrial_projected^* L_F^-1 Rtrial_projected in normalized low coordinates. It is > (1/8-(200/51)10^-26) I, hence >1/9 I. Thus the finite structured block is positive on E85 + F_{85,410}, F=span(e86,...,e410) intersect H0. There are 161 odd and 162 even high directions, giving 408 total finite dimensions (204 per parity), physically orthogonal low/high sectors. W dominates this structured block. This is not E410; the two pole bridges are absent from F. The standalone nonpositive index upper bound 66 is unchanged.

## Exact finite-band rejection for the new trial

Extend ordinary residual U-Lraw X to 600, then project the whole checked band onto its genuine pole-null subspace. The resulting G_600 is a lower bound on the complete residual Gram R^*R, including no assumption about the missing pole bridges. K-(200/51)G_600 has frozen rational negative witnesses: odd <-2/5, even <-1, with upper endpoints -0.42774289703918256 and -1.002702252867781. Thus the particular enlarged trial and scalar penalty fail, regardless of the omitted tail bound. This does not falsify K-R^*L^-1R or original Weil positivity. The full infinite residual was not evaluated after this rejection.

The standard-library independent checker pins both saved result hashes, reruns all 85 shifted integer pivots, directly recomputes negative witness upper endpoints from separate Knew and Gfull intervals (not only their precomputed difference), verifies the interval subtraction and trial Gram upper. No reconstruction of all analytic inputs or numerical solve is claimed.

## Chebyshev residual correction (global analytic theorem; not numerically applied)

Let R=B_V-LX, D=(gamma I-L)/delta=-Q85 P Q85/delta, delta=21/8. Self-adjointness gives spectrum(D) in [-1,1]. For Chebyshev T_m let r_m(t)=T_m((gamma-t)/delta)/T_m(gamma/delta), where gamma/delta=192/175. r_m(0)=1, so q_m(t)=(1-r_m(t))/t is a polynomial of degree m-1. Define X_m=X+q_m(L)R. Its complete residual is R_m=r_m(L)R and R_m^*R_m <= eta_m^2 R^*R, eta_m=1/T_m(192/175). No finite truncation, bridge omission or spectral sampling is used in this analytical bound.

The exact inverse identity gives H_V-B_V^*L^-1B_V = K_m-R_m^*L^-1 R_m, where K_m=K-R^*(2q_m(L)-L q_m(L)^2)R. Hence a sufficient whole-collar condition is K_m-(200/51) eta_m^2 R^*R PSD. Equivalently use the actual R_m Gram, potentially sharper. The scalar penalty applies to a globally contracted residual; the polynomial term retains L-dependent energy. Exact rational recurrence T_{m+1}=2*(192/175)T_m-T_{m-1} proves squared contraction eta_8^2<1/250, eta_12^2<1/8000. Degree 4,8,12,16 penalty coefficients are approximately 0.44704135, 0.01432683, 0.000434144, 0.0000131332 respectively. This does NOT guarantee K_m stays positive: actual polynomial energy and residual moments are still required. It is an application of standard Chebyshev functional calculus, not a novel general polynomial theorem.

## Reproduction

Extract archive preserving cr_alp_045, cr_alp_046, cr_alp_049, cr_alp_050 layout. Run `python cr_alp_050/check_saved.py` (standard library) for saved matrix replay. Run `OPENBLAS_NUM_THREADS=1 python cr_alp_050/audit.py` (numpy, mpmath required) for entry and candidate generation; the numerical solve may vary slightly by BLAS, so regenerated result hashes may differ. The replay checker pins the archived frozen inputs, not arbitrary regenerated candidates. The archived log is the executed original candidate run; predecessor pin checks were added afterward. The final infinite polynomial Schur comparison, fixed-collar positivity, arbitrary support, Lean and RH remain OPEN. Next: calculate rigorously bounded matrix-valued polynomial moments on the original projected translation operator, first test them against both new exact finite witnesses.
