# CR-ALP-123R — 123O obstruction and independent second correction specification

123O successfully completed run 45925861807, job 140266876449, durable commit d4846fc, result blob 1928d1ed6994173e21497660e987bc86e5335515, SHA256 869a0d0a8d8084c396a4e0d2cb956ed404df83d0d04fe0ee51baa380c0698019. Direct whole JSON retrieval was unavailable (empty/truncated connector payload); the following diagnostics use successful runner logs, not independent replay.

E365:
- cosh: E=2.3717871627184496e-38, U=5.716126886186076e-38, tail upper=2.7803906759745824e-19, scalar matrix-order loss/E=558302.8403049805.
- sinh: E=3.761150488339508e-33, U=1.831563576248628e-32, tail upper=1.694018835374374e-16, loss/E=53309.62326382556.

Thus U/E is about 2.41 and 4.87. Even removing every finite coupling penalty cannot make the currently emitted upper-charge sufficient receiver pass. This is not a physical inverse-response lower bound, and does not disprove positivity. 123P activated at 88088746c37ba0550b40773df43e5a19881d65f7 to diagnose the actual matrix directions, not expected to promote positivity with those U inputs.

## Proposed exact correction receiver, preserving the 123N2 source trial

Set q=-y, h=(I-S85)b, Qe_j=e_j-c_j h/st. On a high index bank through cutoff J1 (suggest 850, within the historical correction support), form the real SPD compressed C matrix:
 
 C_ij = P_ij - c_i/st * <h,P e_j> - c_j/st * <h,P e_i>
        + c_i c_j/st^2 * <h,P h>
        + delta_ij d_i - c_i c_j/st*(d_i+d_j)
        + c_i c_j/st^2 * <h,Dcap h>.

This requires only prime columns and exact capped diagonal Parseval, not dense Gamma coefficients among high modes. Choose an approximate numerical solution Cbank w=e_bank, then round to exact rational coefficients. Never accept the numerical solve as a proof. Compute the new q'=q-w and tau' exactly from its high precursor. Independently recompute e' via the centered complete coefficient receiver using f unchanged and q' in place of q; retain complete pole projection in all terms. At fixed f the new joint energy is enclosed using

 E'=E-2 sum w_j e_j+sum w_i C_ij w_j,

where e_j and C_ij are actual outward interval enclosures, and the w_j are exact rationals. If E' has negative upper endpoint this is a lower-comparator trial only, not an original Weil witness.

The conditioned finite covector is no longer negligible: for each original E85 vector v_i, R_i(w)=sum_j w_j (W(v_i,e_j)-c_j/st W(v_i,h)). Apply the frozen V conditioning to this signed vector. The new finite residual is old conditioned residual minus V*R(w); preserve the old certified norm as a small uncertainty and the new signed midpoint as a vector. Then b_(y') = old finite residual - V*R(w) - X*e'. Apply J permutation before computing exact inverse quadratic using K,N. All signs matter; do not sum separate norms of V*R(w) and X*e' when signed cancellation is available.

The finite-bank C solve does not require the nearly-null low pole optimization. Increasing J1 means the geometric tail remainder order also needs to increase, with explicit rho=(max nu_j / (15*J/11))^2<1. Reuse centered denominator expansion, with e' orthogonal to h exactly. Evaluate the weighted pole subtraction using the new e'. Store prefix, tail, C solve defect, signed finite response, exact trial coordinates and updated energy. Compare against the same frozen 123O trial, independently evaluate all charges. If E' is not positive or the bound is inconclusive, no sign promotion.

## Scope

This is a concrete receiving specification, not an executed second correction. 123P is the sole submitted additional workflow. No new local certificate or global support claim.
