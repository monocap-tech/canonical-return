# CR-INT-T1D-S — Actual parity Ritz lower witnesses and the logarithmic Galerkin-tail barrier

**Canonical Return · 2026-11-23 (America/Los_Angeles).**  
**Disposition:** **DERIVED / ANALYTIC, finite exact-rational and symbolic AUDIT-PASS**, based on the original signed Weil smooth-core results T1C, T1D-F/K/R. This supplies **genuine one-sided lower spectral bounds** for the *actual* physical-complement Birman–Schwinger operator and a quantitative necessity theorem for the **particular** crude R.11 Galerkin tail estimator. It **does not** compute the leading low-frequency eigenvalue, prove a spectral upper bound below 1, establish complete complement positivity, identify Connes/Sonin domains, certify Lean, or prove RH. Draft PR #2 remains under review.

## 1. The exact spectral question, unchanged

Use the original signed Weil form of Suzuki, [arXiv:2606.09096v4, §1.1 and §2.4](https://arxiv.org/html/2606.09096v4), on the compact smooth, doubly Mellin-pole-null source
\[
\mathcal D^0=\{f\in C_c^\infty(-2,2):\int f(x)e^{\pm x/2}dx=0\}.
\]
Preserve the physical \(L^2(-2,2)\) inner product (linear first argument), the fixed T1C width \(\ell=10^{-4}\), and the exact two-dimensional selected space \(E=\operatorname{span}\{u_+,u_-\}\). Let \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\) on the window, \(\mathscr H_F=(E\oplus P)^\perp\) and \(F_s=\mathcal D^0\cap E^\perp\). The rank-two source Weil cross penalty is \(\mathcal R=\sum_\sigma |r_\sigma\rangle\langle r_\sigma|/\mu_\sigma\ge0\), from T1D-E/F. The genuine prime operator \(A\) is positive-*kernel*, self-adjoint, and enters the form as \(-A\), without changing sign.

As established analytically in [T1D-R](CR-INT-T1D-R.md), the positive form \(d[y]=\gamma[y]+6\|y\|_2^2\) defines \(D_F\ge5/8\) with compact resolvent. With \(C_F=P_FA P_F+6I+\mathcal R\), the compact self-adjoint operator
\[
\mathscr B_F=D_F^{-1/2}C_FD_F^{-1/2}
\]
obeys \(q|_{F_s}-\mathcal R\ge0\iff\sup\sigma(\mathscr B_F)\le1\). It decomposes exactly into parity blocks. **The remaining task is a valid upper enclosure below 1, not another finite positive compression.** Suzuki's underlying finite-window compactness is canonical; our specific Schur correction is an unreviewed source-level derivation.

## 2. Two authentic parity Ritz witnesses, not an invented discretization

Reuse the four genuinely smooth, physically orthonormal, exactly two-pole-null complementary probes from [T1D-K](CR-INT-T1D-K.md): \(v_{h,\sigma}=(y_{h,R}+\epsilon_\sigma y_{h,L})/\sqrt2\), \(h=2,4\), with \(\epsilon_+=1\), \(\epsilon_-=-1\) and \(y_{h,R}=T_{\log h}e_0/\sqrt{J N_y}\), \(y_{h,L}(x)=y_{h,R}(-x)\), \(N_y=22/3\). Their original mollified bump is the fixed T1C bump, not a toy replacement. Define the **unit** same-parity source
\[
\boxed{w_\sigma=\frac{v_{2,\sigma}+v_{4,\sigma}}{\sqrt2}\in F_s,\qquad \|w_\sigma\|_2=1.}\tag{S.1}
\]
It has exactly 12 disjoint compact smooth source cells; all Mellin pole moments vanish. Both selected directions are physically orthogonal to it. It is a finite *trial* vector in the actual infinite-dimensional operator's form domain.

The original prime-power cell Gram from T1D-K, independently checked here by **exhaustive exact rational ratios of all 12 cell centers**, yields
\[
\boxed{\begin{aligned}
\langle Aw_\sigma,w_\sigma\rangle={}&-\frac7{11}\log3
-\epsilon_\sigma\frac{2s}{N_y}\frac{\kappa_4+\kappa_{16}}2
+\kappa_2-\epsilon_\sigma\frac{2s}{N_y}\kappa_8,\\
& s=\frac4{\sqrt3},\quad \kappa_{2^j}=\frac{\log2}{2^{j/2}}.
\end{aligned}}\tag{S.2}
\]
This retains **all** prime-power pairings for this source (\(2,3,4,8,9,16\)); no approximate cancellation or missing non-ternary prime has been assumed. The rational audit obtains strict numerical enclosure diagnostics
\[
\boxed{\begin{aligned}
5.472947968&<6+\langle Aw_+,w_+\rangle<5.472947970,\\
6.109076352&<6+\langle Aw_-,w_-\rangle<6.109076354.
\end{aligned}}\tag{S.3}
\]
(Endpoints above are deliberately wider than the printed checker outward intervals.) Because \(\mathcal R\ge0\), the **full** \(\langle C_Fw_\sigma,w_\sigma\rangle\) is at least the corresponding (S.3) numerator; the Schur penalty is **retained**, not silently set equal to zero.

For the denominator, import the exact T1C source enclosure \(g_0/J<13.44129808\). T1D-J/K bound each translated 3-cell Gamma self correction by \((4s+2)\ell/N_y\), each reflected Gamma correction by \(3\ell(2+s)^2/N_y\), and the parity-\(2\)/parity-\(4\) cross Gamma correction by \(10\ell(2+s)^2/N_y\). Thus
\[
\gamma[w_\sigma]\le \frac{g_0}{J}
+\frac{\ell}{N_y}\big[(4s+2)+13(2+s)^2\big]
<13.44129808+\frac1{100}.
\]
In particular the true shifted Gamma denominator satisfies
\[
\boxed{0<d[w_\sigma]=\gamma[w_\sigma]+6<\frac{39}{2}.}\tag{S.4}
\]
The conservative inequality is shape-aware and applies to the original complete Gamma source, including its diagonal singularity. The exact-rational check uses \(s<7/3\), \(N_y>7\), \(\ell=10^{-4}\), and T1C's source-bound import.

For \(z_\sigma=D_F^{1/2}w_\sigma\ne0\), the **actual** Birman–Schwinger Rayleigh quotient is
\[
\frac{\langle\mathscr B_{F,\sigma}z_\sigma,z_\sigma\rangle}{\|z_\sigma\|_2^2}
=\frac{\langle C_Fw_\sigma,w_\sigma\rangle}{d[w_\sigma]}.
\]
Equations (S.3)–(S.4) and \(\mathcal R\ge0\) imply strict, independently algebraically auditable **one-sided spectral bounds**
\[
\boxed{\sup\sigma(\mathscr B_{F,+})>\frac7{25}=0.28,\qquad
\sup\sigma(\mathscr B_{F,-})>\frac{31}{100}=0.31.}\tag{S.5}
\]
These are **genuine source-defined positive compact Ritz witnesses**, not computed leading eigenvalues. They are **highly localized rather than genuinely low-frequency**. Therefore their bounds are too far below the critical threshold 1 to certify or falsify (R.8). They also do not bound the top of either spectrum **from above**.

## 3. Why the old Galerkin tail bound alone cannot be the practical stopping rule

The important negative-method result is about **the explicit error estimator printed in T1D-R**, *not* about the true spectral approximation error. Conditional on T1D-I and T1D-Q's \(\|C_F\|<15\), R.11 uses
\[
\boxed{\mathcal E_N=30\sqrt{\frac8{5d_{N+1}}}+\frac{15}{d_{N+1}}},\tag{S.6}
\]
where \(d_{N+1}\) is the \((N+1)\)-st eigenvalue of \(D_F\). A numerical use of this stopping estimator with \(\mathcal E_N\le1/2\) **requires** \(d_{N+1}\ge5760\), already from its first nonnegative summand. The actual error could be smaller; this is strictly a limitation of this **particular sufficient bound**.

We can give an entirely source-valid *upper* bound for that gamma eigenvalue. Let \(S_M\) be the span of the first \(M=N+5\) normalized Dirichlet sine eigenfunctions on \((-2,2)\). Its intersection with \(\mathscr H_F\) has dimension at least \(M-4=N+1\), because there are only four independent physical constraints. Although the raw Dirichlet sine functions are not compactly supported smooth tests, their constrained combinations lie in the **form completion** \(\mathscr V_F\): approximate in \(H_0^1\) by compact smooth functions and use four fixed compact smooth moment correctors to preserve all four constraints. The logarithmic Gamma form norm is dominated by \(H^1\), making this approximation valid; the approximation is not asserted to preserve literal compact smoothness at intermediate noncore stages.

Every \(f\in S_M\) satisfies \(\|f'\|_2\le(\pi M/4)\|f\|_2\). Import T1D-F's uniform **whole-line Gamma multiplier upper estimate** \(m_\Gamma(\xi)\le\log(1+|\xi|)+40\). Jensen and Plancherel give
\[
 d[f]/\|f\|_2^2\le 46+\log\bigl(1+\pi M/4\bigr).
\]
The Courant–Fischer min–max principle therefore proves the useful independent inequality
\[
\boxed{d_{N+1}\le46+\log\left(1+\frac{\pi(N+5)}4\right).}\tag{S.7}
\]
Combining (S.6)–(S.7), the **specific R.11 estimator cannot fall to 1/2** unless
\[
\boxed{N\ge\frac4\pi\bigl(e^{5714}-1\bigr)-5.}\tag{S.8}
\]
This is an astronomically large rank, enormously beyond even T1D-Q's millions of conditioned prolate modes. The exact-rational checker verifies the \(5760\) threshold and proves that \(N=10^6,10^8,10^{12}\) cannot possibly satisfy this **particular estimator**. It does not numerically approximate an eigenvalue of \(D_F\) and does not rule out a much sharper *operator-adapted* tail estimate or an efficient finite spectral proof.

## 4. What has and has not changed

- **Established analytically:** the two real parity trial modes (S.1) exist inside the original smooth Mellin-pole-null physical complement; the original prime form yields exactly (S.2); the true Gamma denominator obeys (S.4); the rank-two Schur term is retained by its sign. Therefore the **actual** compact parity spectra have the lower bounds (S.5), though their leading eigenvalues are not identified.
- **Source-level method limitation:** the min–max bound (S.7) makes the previously advertised R.11 **uniform operator-norm tail estimator** infeasible as a practical sub-1/2 stopping rule. This is a provable barrier for **that estimator**, not a no-go theorem for compact spectral certification or Weil positivity.
- **Exact audit scope:** the companion `src/audits/canonical_return_t1d_s_ritz_tail.py` symbolically enumerates all 12 pairwise-distinct source centers and checks their exact prime-power ratios, certifies rational enclosures for both \(6+\langle Aw,w\rangle\), the Gamma correction budget, and the min–max/tail arithmetic. The Gamma source estimates, functional-analytic compactness and form-core approximation are mathematical imports with stated proofs; the checker does not formalize them or calculate the true Ritz top.
- **Still OPEN / NO TRANSFER:** a certified **upper** Ritz enclosure below 1 on either parity block, a realistic and rigorous upper spectral-tail error, complete original pole-direction control, the correct Suzuki/Connes semilocal intertwiner and H2′/R1/R2, full Weil positivity, or RH.

**Next CR-INT-T1D-T:** build a spectral-tail estimate using *actual operator structure*: compress the joint arithmetic-plus-Schur operator relative to the Gamma spectral projector, rather than bounding every cross block by \(\|C_F\|\) and \(d_1\). Quantify \(\|P_{>N}C_FP_{\le N}\|\) or a source-defined commutator \([D_F,C_F]\) on the appropriate form core; unless an extra smoothing or almost-commutation property is established, do **not** advertise computational tractability. In parallel, use smooth **broad-support, lower-frequency** pole-null trials to investigate whether either leading compact eigenvalue approaches the critical level 1. Neither the present concentrated witnesses nor the raw R.11 tail can decide it.