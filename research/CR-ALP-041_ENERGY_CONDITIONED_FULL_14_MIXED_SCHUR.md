# CR-ALP-041 — Energy-conditioned rational congruence closes the fourteen-mode mixed Schur gate

**DRAFT / NONPROMOTING.** New conditional classical analytical deduction with executed exact rational arithmetic. No new Lean proof, complete collar positivity, arbitrary-support extension or RH claim.

## 1. Source and scope

Same original fixed collar a=11/10, two-pole-null physical space H0, signed Weil form W, and compressed original Gamma-minus-seven-prime operator T as CR-ALP-038/039/040. The original source is `Zeta23.EF.literatureRHS`, Weil-Lab `Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`. The seven prime powers are 2,3,4,5,7,8,9 with authentic negative paired translation terms. Both original boundary anchors remain included. Use CR-ALP-027's analytical W(q,q)>=1/4||q||² on Q15000 H0 in the log form domain, CR-ALP-038's whole-tail interior source estimate (38.3), CR-ALP-039's complete anchor bounds 1/25 and 1/19, and CR-ALP-040's original signed lower matrix and physical Gram. These are inherited analytical premises, not machine-proved continuum identities. The present audit checks the numerical implication only. No exact high-sine variance matrix is claimed to have been computed.

Pinned CR-ALP-040 files (SHA256):
- `CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv`: `7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c`
- `CR-ALP-040_initial_rational_intervals_1e25.txt`: `6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9`

## 2. Energy-conditioned Schur criterion

For either parity write p in the ordered original basis (w3,w5,w7,w9,w11,w13,u1) or (w4,w6,w8,w10,w12,w14,u2). Set sigma=10^-14 and K=H_lower-sigma M. Construct real rational coordinate vectors v_i by Gram-Schmidt-like elimination on the exact midpoint of K, rounding coordinates to 10^-18 and scaling each vector by its maximum absolute coefficient. Their triangular shape with nonzero diagonal guarantees invertibility; approximate orthogonality is NOT assumed in the proof.

Let G_ij=v_i^T K v_j. Enclose every entry by exact Fraction arithmetic from the pinned symmetric matrix intervals. If G_ii has lower bound l_i and |G_ij| has upper bound b_ij, use

    d_i = downward_round(l_i - sum_{j != i} b_ij).

Then G >= diag(d_i) by 2|z_i z_j| <= |z_i|²+|z_j|². This also holds for complex coordinates. All fourteen d_i are strictly positive by exact rational checks. For p=sum z_i v_i,

    W(p,p) >= sigma ||p||² + E(z),  E(z)=sum d_i |z_i|².

For a transformed real source v_i = sum x_k w_k + y u_kappa, the interior part is sum d_j^sine e_j, where d_k^sine=x_k and d_kappa^sine=-sum x_k lambda_k. Use the signed moment D=sum d_j^sine nu_j and C_s=sum |d_j^sine| nu_j^s (s=0,2,3). Crucially, combine the coefficients before bounding their absolute values. With a=11/10, S_a<177/50, N=15000, the inherited CR-ALP-038 bound is

    ||Q_N T(interior)|| <= A/sqrt(3N³)+B/sqrt(5N⁵),
    A=(2S_a/a+201*pi/(400a))|D|(2a/3)^2,
    B=(S_a C_0/(2a)+8S_a C_2/a+201 C_3/(200a))(2a/3)^3.

The derivation extends directly to the new indices 13 and 14: the required separation nu_i/nu_j<1/2 holds for all j>N and i<=14. Add |y| epsilon_kappa, epsilon_1=1/25, epsilon_2=1/19, using the original CR-ALP-039 anchor bounds. Enclose r_i>=||Q_N T v_i|| with exact interval D, C_s, pi<22/7 and integer square-root lower bounds for the denominators. Round r_i upwards to 10^-25.

By triangle inequality and weighted Cauchy-Schwarz,

    ||Q_N Tp||² <= (sum r_i |z_i|)^2 <= rho E(z),
    rho=sum r_i²/d_i.

No assumption of orthogonal residual columns is made. Unlike the previous rank-one-plus-diagonal penalty, the new comparison is conditioned in the signed energy metric.

## 3. Executed result

The exact audit proves, for the odd-indexed sine sector, rho<13/125=0.104, and for the even-indexed sine sector, rho<79/1000=0.079. Both are strictly below 1/8. With Young's inequality and the existing positive tail theorem,

    W(p+q,p+q) >= sigma||p||² + (1-8 rho) E(z) + (1/8)||q||².

Since the two reflection parities are orthogonal for the original signed form, and since W(q,q)>=1/4||q||² for pure-tail q, this gives the strict conditional theorem

    W(p+q,p+q) > 10^-14 ||p||² + (1/8)||q||²

for every nonzero p+q with p in E14, q in Q15000 H0 intersect the original log form domain. This resolves the CR-ALP-040 conservative odd-indexed comparison failure. That earlier negative witness remains a valid falsifier of its particular sufficient comparison, not of the original Weil form.

## 4. Reproduction and open obligations

Run `python cr_alp_041_energy_conditioned_schur_audit.py`, then `python cr_alp_041_energy_conditioned_schur_audit.py --verify` in a directory containing the two pinned input interval files. The second invocation rechecks the stored basis and budgets without regenerating the basis. Both invocations passed. The script uses only Python standard-library exact Fraction arithmetic; it does not independently reconstruct the predecessor Gamma bank or formalize the inherited continuum theorems.

The entire intermediate shell E15000 minus E14 (physical orthogonal complement within E15000) still has dimension 14,986. None of its general signed interactions is certified here. The next stage should use energy-conditioned coordinates for a larger shell, while improving the source-column residual bounds if the trace budget exceeds 1/8. The trace criterion is sufficient, not necessary; no negative original witness follows from failure. No receiving acceptance ledger is promoted.
