# CR-ALP-122R — Consolidated conditional original fixed-collar theorem and two-pole completion

Research-only, no promotion to an unconditional source certificate. 122P quantitative result at durable 6e04c83, blob 7fb998716bdb5a8f0c178585f1f8ad00e7d3fe37, was reopened. Independent JS BigInt cross-multiplication verifies the saved finite gap, displacement reserve and physical floor inequalities, without reproducing floating calculations:
[
  {
    "parity": "1",
    "physical_floor": "1e-38",
    "exact_cross_products": "PASS"
  },
  {
    "parity": "2",
    "physical_floor": "1e-34",
    "exact_cross_products": "PASS"
  }
]

## 1. Precise theorem and scope

Let a=11/10, I=(-a,a), H=L2(I,C), with zero extension J. Define F_f(r)=integral Jf(x) exp(irx) dx, Elog=(2pi)^-1 integral log(1+|r|)|F_f|^2, D={f in H:Elog(f)<infinity}, and D0={f in D:F_+(f)=F_-(f)=0}, F_pm=integral_I f(x) exp(+-x/2) dx. W is the continuous sesquilinear extension of the literal original signed EF form on supported Schwartz inputs: full Gamma symbol Re psi(1/4+ir/2)-log pi, both pole evaluations, and the negative paired von Mangoldt shifts for 2,3,4,5,7,8,9. No change of sign or measure.

Theorem (conditional computational receiving result). For all f in D0, W(f,f)>=10^-38 ||f||_H^2, under the physical enclosure and common-operator premises below. In addition W(f,f)>= [10^-38/(17+10^-38)] ||f||_D^2. The latter norm is ||f||_D^2=||f||_H^2+Elog(f). In particular the restricted form is a Hilbert inner product equivalent to the log-domain norm.

## 2. Assumptions and proof dependency boundary

A. The classical digamma partial fraction, resolvent and elementary coefficient identities in 024,046,057,066,068 and 122G hold for the original W. 122A gives classical domain and parity decomposition; 122O pins the independently native-certified IV103 source (staging 018f7e59, gate 45723601649, 11/11 strict named axiom receipts), and IV106-E8 prime-lag dictionary. This does not imply that the numeric banks are Lean-certified.

B. In each parity the 049-conditioned projected sine coordinate map V spans E85, with Q its physical orthogonal complement inside H0; the original mixed block R is bounded, the 056 correction X belongs to Q, and the actual tail W_Q>=C=Q[Dcap+P]Q>=51/200 I. The finite operator bounds and cap include all infinitely many remaining sine indices. Classical derivations: 122A,B,C,G and 066/068; exact arithmetic 122C/D/E/F/H. The pinned original signed matrix, pole-centering tail and source projection must refer to this same W, not just matrices with compatible sizes.

C. Physical proxy approximation, ideal quadrature, midpoint fast-path and outward intervals bound the actual original residual and all used moment/response integrals, under 058 VALIDATION.md and the explicit IEEE binary128/MPFR contract of 122I/J/K/M/Q. Replay and exact endpoint arithmetic establish consistency, but not a mechanically proved hardware model or an independent integral enclosure for every historical output. The actual archived enclosures remain an explicit trusted computational premise. The 122N direct diagonal and affine/block majorant derivations, 122B cap correction and 111 response transport must give simultaneous U_l*B_C U_l<=A_l, with B_C=(R-CX)* C^-1 (R-CX). K must lower-bound the actual corrected energy K_C=V*WV-R*X-X*R+X*CX (V*WV denotes the finite form matrix, not application of an unbounded operator to arbitrary vectors).

D. The frozen 121 exact terminal matrices are lower grids of KTK-K on the same common coordinates, T=sum_l w_l U_l A_l^-1 U_l*, weights nonnegative and sum one. 122L and 122P verify positive shifts delta=1e-5/1e-4, K>=lambda I with lambda=1/10, and rational norm bounds. These are exact certificate computations conditional on B and C, not additional empirical positivity evidence.

Proof. For each parity, the contraction factor B_C^(1/2) U_l A_l^-1/2 yields B_C^(1/2) T B_C^(1/2)<=I. KTK-K>=delta I implies T>=K^-1+delta K^-2 >=(1+delta/kappa)K^-1, kappa>=||K||. Hence B_C<=T^-1 <= kappa/(kappa+delta) K. No inverse of B_C is used. Complete the square in q+Xz+C^-1(R-CX)z for f=Vz+q, obtaining W>=m||z||^2+c||q+Mz||^2, m=delta lambda/(kappa+delta), c=51/200, M=X+C^-1(R-CX). The bound ||C^-1(R-CX)||^2<=kappa/c yields ||M||^2<=2||X||^2+2kappa/c. The physical Pythagorean inequality ||f||^2<= (||V||^2+2||M||^2)||z||^2+2||q+Mz||^2 yields the saved quantitative floors. Reflection orthogonality adds the two parity bounds. For the log norm, combine W>=epsilon||f||^2 and W>=Elog-16||f||^2 from 024 using theta=epsilon/(17+epsilon), epsilon=1e-38. Then W>=(1-theta)epsilon||f||^2+theta(Elog-16||f||^2)=theta||f||_D^2. The upper bound follows from form continuity. Q and the square completion need not preserve smoothness.

## 3. Next exact reduction: restore the two pole directions

Conditional on this theorem, the full fixed-collar problem has at most two remaining nonpositive directions over C, one per reflection parity. This is a codimension bound, not the sign of either direction. More precisely define b_e=cosh(x/2), b_o=sinh(x/2) on I (zero extended), and h_s=b_s/||b_s||^2, s=e,o. They lie in D by the endpoint-jump Fourier bound of 122A. The moment m_s(f)=inner(b_s,f) is linear (first-antilinear convention), and m_s(h_s)=1. Within parity every f has a unique decomposition f=q+t h_s with q in D0 of that parity, t=m_s(f). The original pole term is 2(|m_e|^2-|m_o|^2).

Because W restricted to D0 is an equivalent Hilbert inner product by section 1, the bounded functional q -> W(h_s,q) has a unique Riesz representer z_s in D0 of the same parity, W(z_s,q)=W(h_s,q) for all such q. Equivalently W(q,z_s)=W(q,h_s). Then
W(q+t h_s,q+t h_s)=W(q+t z_s,q+t z_s)+|t|^2 sigma_s,
sigma_s=W(h_s,h_s)-W(z_s,z_s).
By parity the mixed e/o term vanishes. Thus full fixed-collar nonnegativity is equivalent to sigma_e>=0 and sigma_o>=0; strict positivity equivalent to both >0. This is an exact Schur reduction in the complete log domain, not a replacement by the first two sine modes. On D the nonpositive index equals the number of nonpositive sigma_s; a negative sigma gives a negative direction h_s-z_s, which can also be approximated by smooth supported inputs using 122A. Zero signs need care on the smooth core (the completed-space null vector need not be smooth). A positive sigma requires a certified upper bound on W(z_s,z_s), equivalently on sup_q [2 Re W(h_s,q)-W(q,q)]. Finite variational trial values are lower bounds on that supremum and cannot certify positivity. For a strictly negative sigma a finite lower trial bound exceeding W(h_s,h_s) suffices. No numerical sign for either sigma is claimed.

A possible next route is to test the two normalized pole directions directly with the exact existing b-b and sine-pole coefficients of 046, then bound the remaining mixed response by a certified coercive comparison rather than use the very small global epsilon blindly. These are new mathematical obligations beyond the pole-null receiver. Neither a two-pole completion nor positivity on this fixed radius would establish arbitrary support.

## 4. Verification status

No new runner. This is a consolidated classical conditional theorem, not a native Lean proof or unconditional acceptance of original numeric enclosures. The two-pole Riesz/Schur reduction is newly derived here; it relies on the conditional pole-null theorem. Source pins, exact arithmetic receipts and remaining trusted numerical premises are separate. The global Weil and RH ledgers remain unchanged.
