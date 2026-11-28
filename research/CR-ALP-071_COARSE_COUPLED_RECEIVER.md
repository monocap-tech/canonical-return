# CR-ALP-071 — coarse coupled receiver, exact conditional pair PASS

DRAFT / NONPROMOTING. Freeze job 139561588240 of run 45715854023 completed successfully. Artifact 15070700067, SHA256 ebbc51e56e67e7438483928b3570a320b77df06b92a574ad2a0e998fdb4b0ead, was downloaded and independently verified (ZIP CRC, all 11 manifest file hashes). First source and trial reproduced exactly. Second source is coordinate 38 (zero-based 37), independent of first. Its native coefficient log encloses proxy norm squared in [29147727291,29623356494]/10^15; the second omitted-group norm upper is <1.306e-13.

## A cheaper receiving gate

Use the identical capped C=L+A of CR-ALP-070, C>=dmin I, 0<=A<=Acap I. For exact original corrections U_i, define a_i=<U_i,A U_i>, r_i=R_i-LU_i, s_i=r_i-AU_i. The actual pair Schur matrix is S=K_L+U*AU-s*C^-1 s. Contractivity of Q85 and annihilation of the low sine and pole components give a_i<=Acap times the squared l2 norm of the complete finite high-sine precursor, without truncating the projected pole tail.

The exact interval bank gives a_1<3.158747, a_2<0.000839132. Since A C^-1 A<=A, the triangle inequality in the C^-1 metric gives B_2=(sqrt(b_2)+sqrt(a_2))^2 as a sufficient upper bound for <s_2,C^-1 s_2>, where b_2=(sqrt(Ghi)+eps)^2/dmin and eps is an upper bound on the second total L2 proxy error.

Use the 070 first corrected inverse upper B_1 and first Schur margin m_1. With m_2=K_L22_lower-B_2 and c=|K_L12|_upper+sqrt(a_1 a_2)+sqrt(B_1 B_2), the true S satisfies S11>=m_1, S22>=m_2, |S12|<=c. This is an entrywise sufficient test, not a claim of a signed matrix Loewner lower bound with off-diagonal -c. If m_1,m_2>0 and m_1*m_2-c^2>0, S is positive definite.

## Independent exact arithmetic result

With the explicit conservative premise eps<=1/1000 and the logged second proxy norm enclosure, all arithmetic is exact Fraction with outward integer square roots. Rounded display bounds:

- second old inverse energy <=0.000156426
- second corrected inverse energy <=0.001720159
- S11 >=0.061214171350
- S22 >=0.997465874915
- |S12| <=0.074068519
- det S >=0.055572901531

Status PASS_CONDITIONAL_PAIR. The determinant test passed in the local runtime on the downloaded native artifact. The compact repository checker is src/audits/cr_alp_071_coarse_receiver.py; its native replay is not yet run. No cross integrations, new response solve, or discarded cross terms are required for this sufficient test.

## Remaining attachment

The 058 note states an inherited 057 operator proxy error of 5e-9 times the witness norm, plus omitted groups. This is arithmetically far below 1e-3 for the second direction. However, the numerical midpoint/forward-error and quadrature envelope of the expanded-source G log must be verified for this direction; do not infer it solely from CI success, cancellation, small omitted-group error or the original witness. The second correction bank must retain the same original-source Q85 construction. Pending correlated physical jobs are left untouched and may provide tighter independent cross estimates, but are no longer logically required by this coarse sufficient gate if its premises are discharged.

Prior 1048-dimensional positivity, conditional nonpositive-index cap <=66, full collar, formal Lean and RH status unchanged. No ledger promotion.