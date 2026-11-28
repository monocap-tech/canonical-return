# CR-INT-T1D-AX — Full arithmetic-slack variational trial

**Disposition:** local computational certificate with directed integer-grid arithmetic. Conditional on AP/AM/AW's analytical source identities, inherited AQ/AS/AW arrays and AV high-subspace certificates, and AR's ordinary binary64 exact-integer-dot model. The large AV inverse and original AQ scan are not rerun here. No Lean, remote CI, independent mathematical review, or whole fixed-window/global positivity. No original negative witness.

## 1. Recover the full arithmetic slack

Keep the original I=(-2,2) physical L² metric, smooth double Mellin-pole-null core, both parity embeddings of AP, and all 24 original prime powers of A (positive-kernel self-adjoint arithmetic translation operator, entering the signed Weil form with a minus sign). Inherit AM's audited norm ||A||<lambda=399/50. With t=1/2 the unused source slack in AW is S=(3/5)(lambda I−A), positive. Define B=S/t=(6/5)(lambda I−A). Then 0<=B and I+B<=21I (since 1+12lambda/5=20.152<21). The original source comparison, with the even constant still excluded, is

    q0[f] >= t(<(D_tot+B)f,f> − ||T*f||²) + G_sigma[f],

where D_tot=D_H+R_L diag(1−alpha_L)^−1 R_L*, T=[T_L,T_H], and G_sigma is AW's positive full Gamma image correction. Thus the additional arithmetic energy has not been counted twice. The high deficit and surplus arrays and AV certificate are unchanged.

## 2. Joint variational inequality with full B

Define Z=T_L*R_L, W=T_L−R_L Z*, c=1/21. For physical r=T_L f+T_H g, set U=[R_L,B^(1/2)] and D0=diag(diag(1−alpha_L),I). For any trial v=(v_L,v_B), the inverse variational inequality is

    <r,(D_H+U D0^-1 U*)^-1 r> <= <r−Uv,D_H^-1(r−Uv)> + <v,D0 v>.

This identity follows by minimizing the right-hand side over v; it is valid even for singular B, since D0 is strictly positive. Choose v_L=Z* f and v_B=c B^(1/2) W f (zero trial columns for high coefficients g). Then r−Uv=(I−cB)W f+T_H g. No operator square root has to be evaluated numerically. On the LL block use D_H^-1<=I and R_L*R_L+diag(1−alpha_L)<=I to obtain

    LL <= M−ZZ*−2c W*BW + c² W*(B²+B)W <= M−ZZ*−c W*BW,

where M=T_L*T_L and B²+B<=21B. Keep the exact trial mixed blocks:

    X_AX = X−Z X_s−c W*B T_H,
    Y_AX = Y−Z Y_s−c W*B R_H,
    LH = X_AX − Y_AX H^-1 C_H*,
    HH = T_H*D_H^-1 T_H.

Here X,Y,X_s,Y_s,H,C_H have precisely AW's definitions. Dropping only a nonnegative correction from the LL trial preserves the joint block upper order. The same inherited bound ||H^-1 C_H*||<16 applies. This is an upper comparator, not an identity for the exact inverse blocks. The arithmetic correction retains the entire signed operator, rather than a finite-rank SPS minorant.

## 3. Signed arithmetic matrix construction

For positive circular indices j,k of the same parity p, define omega_j=j*pi/8, sigma=(-1)^p, and u_j(x)=cos(omega_j x)/2 (even) or sin(omega_j x)/2 (odd) on I. For s=log n, parity makes the two shifts' integrals equal. The matrix entry of their sum is exactly

    (1/4) integral_{-2}^{2-s} [cos((omega_j−omega_k)x+omega_j s)
                               + sigma cos((omega_j+omega_k)x+omega_j s)] dx.

For j!=k, the difference contribution is (2/pi)/(j−k) times

    sin(pi(j−k)/4 + omega_k s) − sin(−pi(j−k)/4 + omega_j s).

For j=k it is (4−s)cos(omega_j s)/4. The sum contribution is sigma*(2/pi)/(j+k) times

    sin(pi(j+k)/4 − omega_k s) − sin(−pi(j+k)/4 + omega_j s).

This retains physical overlap lengths, zero extension, signs, and every actual prime power. The even constant is not evaluated and remains constrained. Multiply by the original kappa_n=log(p)/sqrt(n), sum, and form J0=(lambda I−A)/4. Since cB=(8/35)J0, the low and mixed corrections are computed as W*J0 W, W*J0 T_H, W*J0 R_H using the signed cross matrices and Z. No separate approximation of B² is required.

`ax_arithmetic.py` generates certified sparse phases from AQ's rational Taylor seed and exact complex multiplication method, with phase grid 10^14 and physical entry grid 10^9. The same seed error <2/D and induction e_k<=8k/D for k<=8192 apply at D=10^14; the checked two-digit product error plus rounding is <4/10^9. Every sine/cosine of pi(j±k)/4 is exactly in {-1,0,1}, so the numerator intervals use radius 8/10^9. Rational log/sqrt weights, directed endpoints for 2/pi, directed divisions, and integer square roots enclose every shifted pair; radix splitting at 10^9 handles large positive sqrt weights without overflow. Matrix midpoints and row/column error sums are placed on the 2^25 grid. All accepted products in `ax_full.py` use AR's checked exact integer dot method, with ordinary binary64 as explicitly assumed in AR; all proposal eigenspaces are untrusted. The weighted cross Gram error is bounded by sqrt(max row error * max column error); multiplication uses norm error propagation plus directed floor Frobenius error. The subspace verification is AW's exact Gram and projected row/column check. The exact phase and overlap derivation above is analytical, not formally verified by the scripts. The low and high source arrays are hash-verified; AM's arithmetic norm and high floor auditor are re-executed by the arithmetic module.

## 4. Executed certificates and custody

Use AW's high block bounds beta_H < 996/1000 even and < 992/1000 odd, rechecked against inherited AV exact logs. Choose beta=997/1000 and verify gamma² < (beta−beta_L)(beta−beta_H), where gamma=||P_L X_AX||+16||P_L Y_AX||. The numerical proposal eliminates 25 even or 20 odd leading directions of X_AX X_AX*+100 Y_AX Y_AX*, then selects eigenvectors of the compressed LL proposal below .95 or .98 respectively. Only the exact interval checks accept the subspace.

| parity | absorbed low dim | beta_L diagnostic | gamma diagnostic | gamma² diagnostic | strict budget diagnostic | remaining constraints |
|---|---:|---:|---:|---:|---:|---:|
| even | 230 | .944048514583747 | .006403206140609621 | .00004100104887914076 | .00005295148541625302 | 719 |
| odd | 234 | .9705581819927882 | .010709063895178977 | .00011468404951102592 | .00013220909003605906 | 724 |

All comparisons are exact fractions in the accepted logs, not rounded decimals. Even full corrected direct run `ax_final_even.log` exited zero. Odd full matrix construction and direct 230-dimensional certificate `ax_full_odd.log` exited zero (25 removed leading directions, .95 cutoff); the saved enclosed odd matrix (SHA-256 in `ax_cache_sha256.txt`) was then independently checked for 234 dimensions (20, .98) by `ax_full_check.py`, log `full_o997_20_.98.log`. The full odd rerun with final default flags timed out during mixed matrix construction; `ax_final_odd.log` is explicitly incomplete, not evidence of success. The final script defaults reproduce the selected parameters; the accepted odd count is supported by the earlier completed matrix construction and cached subspace check. Initial rank-64 SPS trial and excessively coarse product-grid experiments are not relied upon. No new matrix upper count is obtained from incomplete runs.

On the original smooth pole-null core with the even constant excluded and T*f restricted to the joint certified coefficient subspace, the inverse bound yields ||T*f||²<=beta<(D_tot+B)f,f> by Cauchy in the positive metric. Hence q0[f]>=(3/2000)<(D_tot+B)f,f>+G_sigma[f]>=(3/2000)||f||². Codimension <=(949−230)+(958−234)=1443, parity 719/724. AW's earlier 1501 count is reduced by 58, comprising 441 low-frequency constraints (including even constant) and 1002 unchanged high spectral constraints. This is only an upper bound on possible negative directions; the complementary sign and mixed interactions are not certified.

## 5. A small but exact bound on the remaining Gamma image correction

AW's image correction on the pole-null source is

    G_sigma[f] = sum_{k>=1} [2 exp(−8a_k)/(1−sigma exp(−8a_k))] |<f,v_k>|²,
    a_k=2k+1/2, v_k=cosh(a_k x) or sinh(a_k x).

The k=0 term vanishes by the two pole constraints. Both parity norms obey ||v_k||²<=2+exp(4a_k)/(4a_k). Thus the trace of the positive sum for k>=1 is bounded by sum [4 exp(−8a_k)+exp(−4a_k)/(2a_k)]/(1−exp(−8a_k)). Using r=1/22000>exp(−10), d=1/2900>exp(−8), one obtains

    0 <= G_sigma[f] <= [(4r²+r/5)/((1−r²)(1−d))] ||f||² < 10^-5 ||f||².

`ax_gamma_bound.py` verifies the exponential inequalities with finite positive Taylor sums and the last bound exactly (4257200/467705332367). This does not imply the correction is irrelevant for true near-null directions (AF), but it shows it is not by itself a large uniform reservoir of positive energy for the current comparator. It is unconsumed in the accepted AX spectral result.

## 6. Reproduction and next gate — AY

Use Python/NumPy/SciPy in the extracted bundle directory. Verify `sha256sum -c array_sha256.txt`, `sha256sum -c aw_sha256.txt`; run `python av_check.py`, `python ax_gamma_bound.py`. Reproduce even with `OPENBLAS_NUM_THREADS=3 PARITY=0 python ax_full.py`; reproduce odd with `OPENBLAS_NUM_THREADS=3 PARITY=1 AX_DROP=20 AX_THRESH=.98 python ax_full.py` (the complete odd rerun at these settings has not yet been observed to finish). To reproduce the accepted split odd workflow, run `OPENBLAS_NUM_THREADS=3 PARITY=1 AX_DROP=25 AX_THRESH=.95 AX_ALLOW_SMALL=1 CACHE=1 python ax_full.py` and then `OPENBLAS_NUM_THREADS=3 PARITY=1 DROP=20 THRESH=.98 TARGET=997 python ax_full_check.py`. The explicit AX_ALLOW_SMALL=1 flag permits the initial 230-dimensional odd check. The cache is included, so its hash and separate subspace verifier can also be executed immediately. The original AV high certificate may be independently regenerated with the bundled AV scripts; only its log checker was rerun here.

**Next AY:** apply a residual-certified, iterated physical arithmetic-slack trial to the *high* deficit coordinates, not merely low ones. In the joint variational identity choose nonzero high trial columns of B^(1/2) times a physical synthesis operator, so the unchanged 1002 high constraints can potentially be reduced. Direct mixed cross Gram of signed A and actual D_H inverse residuals should replace the crude 16 multiplier where useful. The k>=1 Gamma image can be included as extra positive factors, but the above norm bound cautions against expecting a large uniform index reduction from it alone. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV, and RH remain unresolved.