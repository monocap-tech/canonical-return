# CR-INT-T1D-AD — Constrained concentration, relative graph-norm transfer, and an exact finite mixed-block test

**Disposition:** new analytical results, conditional on the explicit AB/AC/F/AA/Y/Z source estimates. Exact-fraction audit executed locally. No computed spectral sign, certified negative witness, Lean/CI, Phase IV, completed Suzuki–Connes transfer, or RH. The original source and signed prime interactions are unchanged.

## 1. Source and alternative low/high split

Retain AC's physical space H=L²(I), I=(-2,2), eight mutually orthonormal constraint carriers c_i (two poles, two selected vectors and four central/satellite vectors), C=span(c_i), J=P_{C^perp}, Z=C^perp. The smooth core is Z_s=C_c^infty(I) intersect Z. On it the second Schur form is h=q_0-R-T, with the full original multiplier k(xi)=Re psi(1/4+i xi/2)-log pi-2 sum_{p^j<=54} Lambda(p^j)/sqrt(p^j) cos(xi log(p^j)). Both Mellin-pole terms vanish by the source constraints, not by omission. Import |k(xi)|<=65+log(1+|xi|), ||R||<1785, ||T||<1073590. Use the original unitary whole-line Fourier transform of zero extensions; R,T are physical bounded extensions as in AC.

Fix B=2^21, theta=2^-38. Write S=S_B for band-to-interval restriction, T_B=SS*. **Change the order of conditioning and spectral decomposition:** let U=JS and T_Z=UU* on Z. Let L=Ran 1_(theta,1](T_Z), H_hi=Z orthogonal to L. For w in H_hi, ||S*w||²=<T_Z w,w><=theta||w||². AA.10–16 actually use only this leakage inequality, not the original uncompressed prolate orthogonality. Thus h[w]>=d_*||w||², d_*=1233/1750, on its smooth constrained core. By min-max, eigenvalues of T_Z are bounded above by those of T_B; AB's rank bound gives dim L<=2693201. Parity is preserved. This is NOT identification of L with AC's projected old eigenspace J E; it avoids its unknown angle entirely. Finite physical moment constraints are independent on C_c^infty(I) by L² density.

For f in L, h_f=U*(T_Z|_L)^(-1) f satisfies U h_f=f and ||h_f||<=theta^(-1/2)||f||. The latter is a **relative bound in the actual physical low norm**, without an inverse of J|E.

## 2. Relative finite Fourier transfer

Reuse AC's auxiliary B-spline window chi, a=241/120, k=14, C_window=3360, N=1380000. Let S_N be restriction of P_N chi S to I. The AC bounds ||S-S_N||²<=U0, ||partial_x(S-S_N)||²<=U1 with U0/theta<2^-80, U1/theta<2^-36 are imported and checked arithmetically anew. Put g=J S_N h_f and z=(S-S_N)h_f. Then ||z||<t||f|| and ||z'||<d||f|| where t=2^-40, d=2^-18. In particular ||f-g||<t||f||.

For the positive form norm E²=integral (65+log(1+|xi|)) |hat(v^0)|², AC's trace/log lemma gives ||z||_E<10^-10||f||. The eight physical carriers obey ||(1-J)v||_E<=275||v|| (AC.10). Hence ||f-g||_E<eta||f||, eta=2^-31.

The derivative operator on the orthogonal carrier projection has norm less than 53000: the two normalized pole derivatives are <1, two selected derivatives <37000 and four central/satellite derivatives <46. Therefore ||f'|| <= (B+53000)theta^(-1/2)||f|| < 2^41||f||. AC.7–8 with t=1, d=2^41, R_cut=2^82, alpha<2^43, log(1+R_cut)<58 gives ||f||_E<22||f||. Here the derivative is interior; zero extension can jump. Consequently, for nonzero f in L,

- | ||f||²-||g||² | < 2^-38 ||f||²;
- |q_prime[f]-q_prime[g]| < 10^-10 ||f||²;
- |q_Gamma[f]-q_Gamma[g]| and |q_0[f]-q_0[g]| < (44 eta+eta²)||f||² < ||f||²/40000000;
- |h[f]-h[g]| < (44 eta+eta²+(1785+1073590)t(2+t))||f||² < ||f||²/500000.

These are relative in ||f||, unlike AC. All constraints and original source coefficients remain exact. g is in the form closure, not a literal compact smooth source.

## 3. New graph-norm estimate for the complete mixed representer

For v in H¹(I) set t_v=||v||, d_v=||v'|| and alpha=2 sqrt(t_v²/4+2 t_v d_v)+2d_v. AC.7 gives |hat(v^0)(xi)|<=alpha/(sqrt(2pi)|xi|). For r>0 write L_r=log(1+r), W_r=65+L_r. Since for xi>=r, log(1+xi)<=L_r+log(xi/r), integrating the squared weight gives

||k(D)v^0||² <= W_r² t_v² + alpha²/(pi r) (W_r²+2 W_r+2).

Indeed the tail integrals of 1, log u and log² u against u^-2 on [1,infinity) are 1,1,2. This bounds the actual whole-line multiplier image in L², not just its quadratic energy. For z normalized by ||f||, AC's t,d yield alpha<10^-5; r=2^40 yields W_r<94. The rational audit gives ||k(D)z^0||<10^-9||f||.

For normalized sharp pole carriers the same inequality with t=1,d<=1, alpha<=5 and r=2^10 gives ||k(D)c_i^0||<74. For smooth selected carriers the symbol bound and log(1+x)<=x give <37065; for central/satellite carriers <111. Cauchy–Schwarz on the eight coefficients gives ||k(D)((1-J)z)^0||<=53000||z||. Define K_h v=J[(k(D)v^0)|_I - Rv-Tv] on Z intersect H¹(I), so h(v,w)=<K_h v,w> for the source core and its appropriate form extension. J is contractive; no physical window jump is silently removed. Therefore

**||K_h f-K_h g|| < (10^-9+(53000+1785+1073590)2^-40)||f|| < ||f||/950000.**

For any fixed high-sector orthogonal projector P_hi, the same bound applies to P_hi K_h(f-g). In particular the square roots of the corresponding residual quadratic forms differ by <=||f||/950000, a two-sided mixed-transfer estimate. This is NOT an upper bound of 1/950000 on the mixed coupling itself. Its size still has to be evaluated.

## 4. A fully finite alternative split (no exact prolate eigensystem required)

Set V=J S_N, T_fin=VV* on Z, and L_fin=Ran 1_(theta/4,infinity)(T_fin). Its rank is <=2N+1=2760001. Let H_fin=Z orthogonal to L_fin. For w in H_fin,

||S*w||=||U*w|| <= ||V*w||+||U-V|| ||w|| <= (sqrt(theta)/2 + sqrt(U0))||w|| < sqrt(theta)||w||.

Thus the **same** AA high coercivity is valid on the finite-split high sector, with no eigengap assumption. Its proof is conditional on AA's original signed-source inequalities. Every physical vector f in L_fin has a band preimage of norm <=2 theta^(-1/2)||f|| under V. This avoids the old quotient angle; it does not assert that the raw overcomplete trigonometric Gram is well conditioned. It is also not an identification with the exact L from section 1.

Explicit implementation: let Phi contain the restricted Fourier columns of AC.12, F=J Phi, G=F*F (the sinc Gram minus the exact rank-eight moment Gram). Let A map a band function to the Fourier coefficients of P_N chi S h; its kernel is the AB.5 Fourier window divided by sqrt(4 pi a). Set D=A A* (finite positive matrix of band integrals). Then V=F A, T_fin=F D F*, and W=F D^(1/2). Diagonalize the finite positive matrix D^(1/2) G D^(1/2). For an orthonormal eigenvector v_i of eigenvalue lambda_i>theta/4, ell_i=W v_i/sqrt(lambda_i) is an exactly normalized physical low vector. The D-weighted coordinates have inverse norm <=2/sqrt(theta). This is a finite mathematical prescription; actual matrix spectral certification and verified quadrature are not performed here.

## 5. Exact mixed Gram and a source-preserving sufficient finite sign test

Use either low/high split; the finite split is operationally preferable. Let ell_i be an orthonormal physical low basis, H the matrix of h on it (Hermitian), and K columns K_h ell_i. With P_hi the orthogonal high projection, the **actual** mixed Gram is

M=K* P_hi K=K*K-H² >=0.

This is an exact identity on Z: the low projection coefficients of K are H. The mixed penalty in the high Schur minimization is at most d_*^(-1) M by high coercivity and Cauchy–Schwarz. Therefore the finite Loewner certificate

**H-d_*^(-1) M >=0**

is sufficient for h>=0 on the original smooth source, conditional on the preceding AA/Schur hypotheses. It is not necessary, and a negative eigenvalue of this conservative lower matrix is not a counterexample. A strictly negative eigenvalue of H itself does imply a smooth negative source after approximation and completing the earlier finite Schur squares.

For an alternative computable upper enclosure of M define whole-line images U_i=k(D)ell_i^0-(R ell_i+T ell_i)^0, D_R=(<U_j,U_i>) and C0=(<U_j,c_i^0>) with consistent column convention. Low vectors and carriers are mutually orthonormal in the whole line after zero extension. Hence M<=D_R-C0* C0-H²; the difference is precisely the positive Gram of U_i restricted outside I. D_R uses the full k² Fourier integrals and the exact rank-six physical residuals, not a prime-only or periodic Gamma substitute. All H and C0 terms must use the same signed source and same physical metric. The resulting sufficient finite inequality H-d_*^(-1)(D_R-C0* C0-H²)>=0 can be evaluated without enumerating an independent high basis, though its conservatism and computational cost may be prohibitive. No numerical matrix enclosure or mixed Gram value is claimed.

Source custody: the low vectors are H¹(I) functions with zero extensions in H^s(R), s<1/2, and all physical constraints hold. For a smooth w in Z_s, the high component after subtracting its low projection is also H¹(I) with the low/high moment constraints. Approximation in H^s by compact smooth functions and correcting the finitely many independent moments with fixed compact smooth biorthogonal functions gives convergence in logarithmic energy (AC's density argument). Thus high coercivity passes to these high components; completing the low/high squares is valid on the original core. Any negative low vector can likewise be approximated by genuine Z_s tests; the first two positive finite Schur blocks give an original smooth negative q_0 direction. No density condition number for a specific explicit witness is asserted.

## 6. Audit and next gate

Run `python3 src/audits/canonical_return_t1d_ad_constrained_concentration.py`. Exact-fraction checks cover tail and trace constants, relative form and graph-norm errors, and the finite-split leakage margin. They do not check functional analysis, Gamma formula, previous source bounds, spectral values, or matrix quadratures.

**Next AD -> AE:** use the finite constrained split or an independent spectral enclosure to bound the actual low matrix H and mixed residual M (or sharper high inverse in place of d_*^-1). A finite candidate vector with an interval-enclosed negative original q_0 is an alternative decisive falsifier. The new finite test is well-defined and avoids AC's old angle obstruction, but no sign of its lowest eigenvalue is known. Conditional negative-index upper bound remains; zero negative directions, unconditioned complement positivity, all pole carriers, completed Suzuki–Connes transfer, Phase IV and RH remain open. Keep research PR #2 draft.