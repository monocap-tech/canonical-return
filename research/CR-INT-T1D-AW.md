# CR-INT-T1D-AW — Joint low/high deficit absorption with restored low surplus

**Disposition:** local directed-integer computational audit PASS in both parities. Conditional on the original AP/AM source identities, certified inherited AQ/AS arrays and AV high-subspace certificates, and AR's stated ordinary binary64 exact-integer-dot model. AV's large inverse is not rerun here. No Lean, remote CI, independent mathematical review, full fixed-window positivity, or RH. No original negative direction found.

## 1. Complete source and the new low census

Keep I=(-2,2), the physical L² metric, the original two-Mellin-pole-null smooth core, the AP circumference-eight parity embedding E and its normalized cosine/sine restrictions u_j (frequency j*pi/8). The even constant j=0 has separate normalization and remains constrained. Retain AP's hybrid lower diagonal d_j=m_Gamma(j*pi/8)−(2/5)a(j*pi/8)−(3/5)lambda, lambda=399/50, with all 24 original prime powers in a, and reference t=1/2. Reexecute AV's low census with the complete digamma positive series and integral remainder; store downward bounds d_j^lo on the 10^8 grid for every j<1024. The result reproduces 59 even and 60 odd positive baseline coordinates, and identifies 452 positive-index deficit coordinates in each parity, plus the excluded even constant. Define upper deficit b_j=(t−d_j^lo)/t where positive. For the 119 certified positive surplus coordinates restore alpha_j=floor(10^8*(d_j^lo−t)/d_j^lo)/10^8. This is a lower bound on the inverse-surplus coefficient. The inherited high deficit/surplus arrays are SHA256-checked, not rescanned. The inherited AV exact logged spectral bounds are strictly below 996/1000 (even) and 992/1000 (odd), checked as fractions; the actual AV auditors were executed in AV, not in AW.

Let T_L and T_H be physical synthesis operators with columns sqrt(b_j)u_j for low and inherited high deficits. Let R_L and R_H be the analogous surplus synthesis operators with weights sqrt(alpha). Set D_H=I+R_H diag(1−alpha_H)^−1 R_H*, D_tot=D_H+R_L diag(1−alpha_L)^−1 R_L*. On tests with c_0=0 in the even sector, the original lower form is at least t(<D_tot f,f>−||T_L*f||²−||T_H*f||²), plus the original nonnegative Gamma image and arithmetic deficit (Section 4). Parseval is in the actual physical metric. No original signed arithmetic term is replaced by an altered weight.

## 2. Exact joint low/high inverse comparison without another large solve

Write H=diag(1−alpha_H)+R_H*R_H, C_H=T_H*R_H; AU/AV give I/24<=H<=I and D_H^−1=I−R_H H^−1 R_H*. Define

    M=T_L*T_L, X=T_L*T_H, Y=T_L*R_H,
    Z=T_L*R_L, X_s=R_L*T_H, Y_s=R_L*R_H.

Use the variational inverse upper bound for adding R_L, with H_L=diag(1−alpha_L)+R_L*D_H^−1 R_L and B0=T*D_H^−1 T, T=(T_L,T_H). For any trial V, writing Z_true=R_L*D_H^−1 T gives

    T*D_tot^−1 T <= B0−V* Z_true−Z_true* V+V* H_L V.

Choose V=[Z*,0], where the zero columns correspond to high deficits. Since H_L0=diag(1−alpha_L)+R_L*R_L<=I by the circular Parseval frame, the trial's LL block is exactly

    M−2ZZ*+Z H_L0 Z*−(Y−Z Y_s) H^−1 (Y−Z Y_s)*
    <= M−ZZ*.

The LH block is exactly

    X−Z X_s−(Y−Z Y_s) H^−1 C_H*.

The HH block is unchanged from AV: T_H*D_H^−1 T_H. This is a joint operator upper bound before restricting either coefficient block. It is not an assertion that low and high original functions decouple. The inherited high weights b_H<10, Parseval, and H>=I/24 give

    ||H^−1 C_H*||² <= 24 ||C_H H^−1 C_H*|| <= 24 ||T_H*T_H|| < 240 < 16².

The second inequality follows from positivity of T_H*D_H^−1 T_H=M_H−C_H H^−1 C_H*. Hence on any low coefficient subspace S_L, the mixed norm is bounded by

    gamma = ||P_L (X−Z X_s)|| + 16 ||P_L (Y−Z Y_s)||.

If the low compression of M−ZZ* is bounded above by beta_L, and the AV high compression by beta_H, the joint upper is controlled by the scalar 2x2 matrix with diagonal beta_L,beta_H and off-diagonal gamma. To certify a target beta, require beta>beta_L,beta_H and gamma²<(beta−beta_L)(beta−beta_H). This does not need the numerical AV high basis or a new solve. An earlier AW run without restored low surplus used Z=0 and certified 198/194 additional low dimensions at target 999/1000.

## 3. Executed fixed-grid spectral certificate

For the final certificate use the original exact signed physical overlap formula for positive indices, directed sqrt/product enclosures from AR/AS, S=2^25 matrix grid, U=2^23 proposal grid. Compute the original low Gram and signed cross matrices on both deficit and surplus sets. Compute the corrected matrices M−ZZ*, X−Z X_s and Y−Z Y_s by exact integer dots and floor rounding, with independently bounded operator errors. In particular ||Z||<8, ||X_s||<4, ||Y_s||<=1 by Parseval and the checked weight bounds; product norm errors and floor Frobenius errors are added before acceptance. Proposal only: remove 15 leading directions of the positive matrix X_corr X_corr*+100 Y_corr Y_corr*, then choose eigenvectors of the compressed low Gram with eigenvalues <=0.98. Neither proposal eigenvalues nor orthogonality are trusted.

For the rounded proposal Q/U verify its Gram error delta<1 with exact integer dot products. Bound its low compressed quadratic matrix by the induced norm (including two-product rounding and the original matrix enclosure) divided by 1−delta. Bound projected cross operators by the checked row/column norm of Q* times each grid matrix, plus enclosure error multiplied by 1+delta, then divided by 1−delta (conservatively bounding the normalization by 1/sqrt(1−delta)). Every exact integer dot uses the inherited AR ordinary IEEE binary64 model with the independently checked Cauchy bound on absolute sums of products <2^53. The row/column sums and scale products have the corresponding checked int64 bounds. No eigenproposal is an acceptance premise.

Final executed log `joint_final.log` (all displayed decimals below diagnostics only, acceptance uses exact fractions):

| parity | new absorbed low dimensions | low beta | mixed gamma | gamma² | strict mixed budget | remaining constraints |
|---|---:|---:|---:|---:|---:|---:|
| even | 201 | .9681687726 | .0050327731 | .0000253288052 | .0000288312274 | 748 |
| odd | 205 | .9653021376 | .0069286379 | .0000480060233 | .0001584893122 | 753 |

Both exact assertions pass with beta=997/1000, beta_H=996/1000 even and 992/1000 odd. Hence the full retained low/high inverse comparator is <=997/1000 on S_L direct-sum S_H. Impose T*f in that subspace and c_0=0 in the even case. The constraints are linear on the original smooth pole-null source. By Cauchy in the D_tot metric, ||T*f||²<=beta<D_tot f,f>. Consequently q0[f]>=(3/2000)<D_tot f,f>>=(3/2000)||f||², with the original nonnegative corrections still available. The conservative count is (949−201)+(958−205)=1501, including the one constant exclusion, 498 other low-coordinate constraints and the unchanged 1002 high-coordinate constraints. Conditional negative index <=1501 (748 even, 753 odd), 406 below AV. This is not an actual negative eigenvalue count, and does not imply positivity on the remaining sector or its mixed interactions. The AV and AW counts are not double-counted.

## 4. Keep the positive source deficits explicit

The complete AP image formula gives an exact positive Gamma correction on parity sigma=+1 or −1:

    G_sigma[f] = sum_(k>=0) [2 exp(−8 a_k)/(1−sigma exp(−8 a_k))] |<f,v_(k,sigma)>|²,
    a_k=2k+1/2, v_(k,+)=cosh(a_k x), v_(k,−)=sinh(a_k x).

This follows by multiplying AP's resolvent difference by 2a_k in each digamma summand. The k=0 moment vanishes on the original two-pole-null source. The full original equality is the hybrid circle form plus G_sigma plus S[f], where S=(3/5)(lambda I−A)>=0 on the physical source. The gamma sum converges on the compact smooth tests by the original resolvent identity. The numerical AW count does not consume either source correction; indeed on the certified complement it gives q0[f]>=(3/2000)<D_tot f,f>+G_sigma[f]+S[f]. No claim of numerical improvement from either deficit is made.

**Next gate AX:** retain positive source energy as finite-rank factors in the same inverse comparison. For example, any physical orthogonal finite-rank P gives S>=S P S/(6lambda/5), since 0<=S<=6lambda/5 I. This retains signed translations through the actual representers S v, with v an orthonormal physical basis of P, and allows additional exact Gamma image moments to be included directly. Compute their physical mixed Gram and use the joint variational formula, rather than enlarging the same low trial block one mode at a time. This is an analytic construction, not a passed source-deficit spectral certificate. Full fixed-window/global positivity, completed Suzuki–Connes transfer, Phase IV and RH remain open.

## 5. Reproduction

In the bundle directory with Python/NumPy/SciPy:

    sha256sum -c array_sha256.txt
    OPENBLAS_NUM_THREADS=1 python aw_low.py
    python av_check.py
    OPENBLAS_NUM_THREADS=3 python aw_joint.py

`av_check.py` only checks the inherited execution logs; to independently regenerate the high certificate run the included AV inverse auditor for each parity, as described in the AV note. The inherited arrays can be regenerated from AQ/AS scans. The new low and joint audits actually executed successfully; no remote CI or Lean claim.