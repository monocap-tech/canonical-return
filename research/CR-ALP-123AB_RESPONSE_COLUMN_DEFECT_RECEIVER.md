# CR-ALP-123AB — Response-column defect receiver for the infinite auxiliary error

123Z/AA run 45931136929 succeeded, durable 02c5ea4. Source joint-repaired.json SHA256 7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897; joint-result.json SHA256 e03296b5a0717c1447adc3bccc48991add2bd4bb880507b6a216490ad50c0145. All 12 sparse asymptotic tests passed. Joint E850 finite penalty/energy 1.0058326421247526e-27 (cosh), 1.5065204312306415e-29 (sinh); U/E 0.4259758159944923 and 0.3777643993925838. Total sufficient loss/E 16568.26915556367 and 613.1034819911827. These are upper-charge diagnostics, not physical response lower bounds. The observed values are from runner logs, not a new independent physical evaluation.

## 1. Replace universal inverse amplification by a verified response-column defect

In 123F notation C is the full bounded comparator on Q85, R is the actual mixed map in original V-conditioned coordinates, S_C=H-R* C^-1 R>=S0, and S0^-1=D=K^-1+N^-1 in the 082 permutation. Let M=C^-1 R. For the normalized joint trial f+q, put y=-q, e=r-Cy, s=a-R*y, E=W(f,f)-2 Re<r,y>+<y,Cy>. Thus a-M*r=s-M*e. This is an exact identity, not an inequality.

Choose any finite bank response Mhat with columns in Q85; e.g. use projected sine columns through E850 and solve C_B Mhat_B=R_B, using the same bounded block factorization as 123W. Define Z=R-C Mhat. Then M=Mhat+C^-1 Z and

    a-M*r = s-Mhat*e - Z*C^-1 e.

Write b=s-Mhat*e. If U>=<e,C^-1 e> and eta>=||D^(1/2) Z*C^-1 Z D^(1/2)||, the complete pole Schur scalar satisfies

    sigma >= E-U-(sqrt(b*D b)+sqrt(eta U))^2.

Proof: apply Cauchy-Schwarz to D^(1/2) Z* C^-1/2 acting on C^-1/2 e, and use S_C^-1<=D. This is precisely the 123F inverse receiving identity, but Mhat is an approximate full mixed inverse, rather than the historical correction X. Taking Mhat=X recovers 123K with the old defect. If the newly selected columns are confined to E850, Mhat*e uses exactly the near-zero E850 projected residual coefficients; no arbitrary tail extrapolation or additional large amplification of that term is required.

## 2. Computable upper bound on eta, including infinite tail

For each column z_i of Z, let u_i be an upper bound on <z_i,C^-1 z_i>. Then |<z_i,C^-1 z_j>| <= sqrt(u_i u_j). Because the transformed Gram is PSD, its norm is bounded by its trace. Therefore the exact rational sufficient certificate is

    eta = sum_ij |D_ij| sqrt_up(u_i) sqrt_up(u_j).

This is an entrywise absolute trace bound, valid without a commutativity assumption. Sharper trace receiver: compute exact signed coefficient vectors z_j on the explicit high prefix in J coordinates, then evaluate sum_j z_j* D z_j / d_j with outward error bounds, where d_j=Dcap_j-21/8. For the remaining tail, use per-column complete L2 norm bounds tau_i and sum_ij |D_ij| tau_i tau_j/dcap. Indeed C>=Q diag(d_j) Q and dropping the pole constraint only increases the inverse bound on its Q-projected columns. Thus

    eta <= sum_(85<j<=2048) z_j*D z_j/d_j
           + sum_ij |D_ij| tau_i tau_j/dcap.

Use upper interval enclosures on the RHS. This notation means the actual eta is bounded by the displayed quantity; the emitted eta certificate is that quantity or a safe upward bound. The signed prefix expression must retain off-diagonal cancellations. An equivalent computation applies a rational LDL decomposition of D to the coefficient lists and accumulates positive weighted squares. No inverse of the actual infinite operator is evaluated.

A further optional receiver uses the actual bank Schur L_B=H-R_B* C_B^-1 R_B or the energy of Mhat, and the defect Gram; S_C=H-R*Mhat-Mhat*R+Mhat*C Mhat-Z*C^-1 Z. If the latter lower bound improves S0, the same exact matrix-order computation applies. It must be audited separately.

## 3. Required size of the new coupling certificate

Even with zero b, acceptance needs eta<E/U-1. At the current joint E850 U/E diagnostics the limiting values are about 1.34755 and 1.64715. With a nonzero b, use the actual exact inequality, not rounded decimal ratios. The previous universal amplification factors were 38893.85865971763 and 1621.9784568821358. This is a demanding target: the new defect columns must actually be small in the D-weighted geometry. The derivation does not guarantee that they are. If eta remains too large, refine the response-column bank independently of the pole trial or tighten U. Neither the existing U nor finite block stationarity proves the sign.

## 4. Custody and implementation gate

Required inputs: exact rational projected response columns Mhat_B in original V coordinates, exact or outward interval R_B and C_B, bank residuals, complete Z coefficient prefix, signed tail moments and pole-column constants for each Z column, original frozen 123W pole trial and e_bank coefficients, actual signed b, K,N and common J permutation. Use original 046 coefficient formulas for R and capped C, as in 123W; never replace R by a certified lower matrix. 123W currently preserves only the aggregate bank residual bound, not signed e_bank coefficients. Its b can conservatively use ||Mhat_B|| and that aggregate bound, or regenerate exact signed e_bank from the frozen rational trial; do not assume stationarity is exact. For Mhat*e, its precursor coefficients pair with actual e_j because e is orthogonal to h.

Next executable gate is response column generation and complete defect enclosure. No new response column has yet been computed, and no new runner activated. Conditional pole-null index cap 0, unrestricted collar index cap 2; no global positivity or RH claim.
