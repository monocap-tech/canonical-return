# CR-CC-001 — Connes–Consani: precise Weil-window carrier, Sonin metric distortion, and pole-null phase probe

**Canonical Return / independent investigation / 2026-11-22 America/Los_Angeles (2026-11-23 UTC).** Classification: **RECOVERY** of native arithmetic form, **EXACT ANALYTIC SOURCE-ADMISSIBILITY TEST**, **QUANTITATIVE FINITE-PLACE METRIC BOUNDS**, **CONDITIONAL SCHUR RE-ENTRY CANDIDATE**. This is **not** a Lean kernel proof, RH proof, or new unconditional global positivity theorem. The originals at monocap-tech/weil, weil-lab, provenance and hodge-rh are unchanged.

## 1 — Primary papers and scope

- Alain **Connes and Caterina Consani**, *Weil positivity and Trace formula, the archimedean place*, arXiv:[2006.13771v1](https://arxiv.org/html/2006.13771), Theorems 1 and 11: in the prime-free archimedean window, a Weil inequality is compared against a nonnegative Sonin-compression trace. The authors require support and Fourier/pole vanishing conditions; this does **not** prove global Weil positivity.

- Alain Connes and Caterina Consani, *Quasi-inner functions and local factors*, arXiv:[2008.10974v1](https://arxiv.org/html/2008.10974), Proposition 5.24: Sonin space becomes the kernel of the diagonal block \((\rho_\infty)_{22}\). Finite place-set products of local-factor ratios are quasi-inner; their corresponding Sonin spaces form a filtered family.

- Alain Connes, Caterina Consani, **Henri Moscovici**, *Zeta zeros and prolate wave operators: Semilocal adelic operators*, arXiv:[2310.18423v2](https://arxiv.org/html/2310.18423), §4.6–4.7 eq. (57) and Theorem 4.6: an explicit Hilbertian isomorphism \(\theta_S\) between archimedean and semilocal Sonin spaces.

- **Newer primary source**, Ronnie Andrews, Alain Connes, Antonio La Gatta, *The uncertainty operator, the Weil ground state and the zeros of zeta: a numerical study and conditional positivity of the Weil quadratic form*, arXiv:[2611.13129v1](https://arxiv.org/html/2611.13129), posted **November 9, 2026**, §2.1 eqs. (6)–(8), §8.2 Theorem 8.1, §8.5 Theorem 8.14 (**conditional**), §8.6 Lemma 8.16, Proposition 8.18 (**conditional**). This 2026 paper is vital to distinguish current Connes-related results from the September 2026 original RH reference ledger.

- Weil-Lab [IV030 source-pinned phase witness](https://github.com/monocap-tech/weil-lab/blob/1061a331c195e30e85dd8eb35a6eff1c794cf054/docs/IV030_PHASE_IV_ORIGINAL_PRIME_PHASE_FIDELITY_NO_GO.md); original RH [PRED-SEAM-17 exact kernel typing](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_17_SONIN_EQUALITY_20260914.json).

## 2 — Exact COMMON SOURCE CARRIER: the logarithmic isometry

Write \(\lambda=e^a\). The map \(L_a:C_c^\infty((-a,a))\to C_c^\infty((\lambda^{-1},\lambda))\) given by \((L_a\phi)(u)=\phi(\log u)\) preserves the physical Hilbert metric:
\[\boxed{\int_{\lambda^{-1}}^\lambda |L_a\phi(u)|^2\,\frac{du}{u}=\int_{-a}^a|\phi(x)|^2dx.}\]
Moreover, for \(f=L_a\phi\), \(g=L_a\psi\), and multiplicative involution \(g^*(u)=\overline{g(u^{-1})}\),
\[\boxed{(f*g^*)(e^t)=\int_{\mathbb R}\phi(x)\overline{\psi(x-t)}dx=C_{\rm lab}(\psi,\phi;t).}\]
This identity fixes the source's **first-conjugate** Weil-Lab convention and the **linear-first** Connes pairing. The 2026 paper's real-symmetric \(F(x)=\int g(u)h(ux)d^*u\) reduces to the same correlation after argument interchange when both inputs are real. The involution is not complex-conjugate linear in the first original lab slot; slot order must be retained.

## 3 — Full scalar signed Weil form matches, but not an unrestricted operator intertwiner

Under \(u=e^t\), the Connes–Andrews–La Gatta 2026 §2.1 decomposition
\[\Psi(F)=W_{0,2}(F)-W_{\mathbb R}(F)-\sum_p W_p(F)\]
becomes **exactly** the original Zeta23 / Weil-Lab scalar formula on smooth correlations: \(W_{0,2}\) supplies the two exponential pole moments, \(\sum_pW_p\) supplies \(\Lambda(n)n^{-1/2}[k(\log n)+k(-\log n)]\), and \(-W_{\mathbb R}\) supplies the complete Gamma/digamma term. This rests on [CR-SZ-002's checked digamma normalization](https://chatgpt.com/space/page_8c63ddb9ebb481919002432bc72dab3b) and eqs. (6)–(8) of the 2026 Connes paper. No further independent spectral equivalence is implied.

Support \(f\subset[\lambda^{-1},\lambda]\) implies \(f*f^*\) is supported within \([\lambda^{-2},\lambda^2]\): **only prime powers \(n\le\lambda^2=e^{2a}\)** can contribute. Thus n=9 becomes possible exactly for \(a>\log3\) for interior smooth supports, agreeing with [CR-SZ-003](https://chatgpt.com/space/page_f22d2c1da1f0819181937662f715a667).

**Strongest lawful conclusion:** every finite Gram compression of the ORIGINAL smooth-core signed Weil form has the same entries and physical \(L^2\) metric as the corresponding multiplicative-window compression, provided the test objects and slot conventions match. This is an exact form-source identity. The semilocal Sonin representation and completed operator domains require additional maps.

## 4 — The trace-formula pole gate rules out direct use of the old IV030 two-bump span

Theorem 8.1 of the 2026 Connes paper is stated for **real** smooth \(f\) and gives a clean trace representation of the Weil form after both Mellin pole moments vanish:
\[\mathcal M f(1/2)=\mathcal M f(-1/2)=0.\]
Under the logarithmic map these are exactly \(M_\pm(\phi)=\int\phi(x)e^{\pm x/2}dx=0\). Take the original IV030 positive bump \(\phi\) and its \(\log9\)-translate \(\psi\). Put \(A=M_+(\phi)>0\), \(B=M_-(\phi)>0\). Since \(\log9=2\log3\),
\[\begin{pmatrix}M_+(\alpha\phi+\beta\psi)\\ M_-(\alpha\phi+\beta\psi)\end{pmatrix}=\begin{pmatrix}A&3A\\B&B/3\end{pmatrix}\begin{pmatrix}\alpha\\\beta\end{pmatrix}.\]
Its determinant is \(-\tfrac83AB\ne0\). **Exact obstruction:** the two-dimensional IV030 span intersects the *two-pole-null test subspace* only at zero. Its (1, −1, i) phase witnesses cannot be inserted directly into the no-pole trace identity of Theorem 8.1; using it without pole correction would be invalid. This does not prevent computing their full Weil form by the explicit formula, or using a trace identity that separately retains the pole term.

## 5 — A new explicit FOUR-CELL pole-null phase calibration

Construct \(4\) disjoint translated copies \(\phi_j(x)=\phi(x-jd)\), \(j=0,1,2,3\), of a single real nonzero \(C_c^\infty\) bump, with \(d=\log3\); center them symmetrically near \((-3d/2,-d/2,d/2,3d/2)\). Require **\(a>3d/2\)** and choose the bump sufficiently narrow to exclude every cross-prime-power shift except \(d,2d,3d\), i.e. n=3,9,27. Put \(r=\sqrt3\), \(s=r+r^{-1}=4/\sqrt3\), and for \(\beta\in\{1,-1,i\}\) define the coefficient polynomial
\[p_\beta(z)=(z-r)(z-r^{-1})(1+\beta z)=1+(\beta-s)z+(1-s\beta)z^2+\beta z^3.\]
Let \(v_\beta=\phi_0+(\beta-s)\phi_1+(1-s\beta)\phi_2+\beta\phi_3\), \(J=\|\phi\|_2^2>0\). For either sign, translation multiplies the pole moment by \(r^{\pm j}\); hence \(M_\pm(v_\beta)=M_\pm(\phi_0)p_\beta(r^{\pm1})=0\) **exactly**, for all three complex phases.

Every coefficient is nonzero for the chosen \(\beta\), so all three tests have exactly the **same 4-component physical support** and the same arithmetic support-interaction graph. The squared \(L^2\)-norm is \(\|v_\beta\|^2=J N(\beta)\) with
\[N(\beta)=\frac{44-16\sqrt3\,\operatorname{Re}\beta}{3}>0.\]
Thus \(h_\beta=v_\beta/\sqrt{J N(\beta)}\) has **unit mass**, shared support, shared graph and **both pole moments zero**. The original *positive-named* von Mangoldt correlation sum is
\[\boxed{P_{\rm original}(h_\beta,h_\beta)=\log3\,\frac{-7+3\sqrt3\,\operatorname{Re}\beta}{11-4\sqrt3\,\operatorname{Re}\beta}.}\]
Its three exact values divided by \(\log3\), ordered \(\beta=1,-1,i\), are
\[\frac{-41+5\sqrt3}{73},\qquad \frac{-41-5\sqrt3}{73},\qquad -\frac7{11},\]
which are **pairwise distinct**. The genuine signed Weil prime contribution is their **negative**. The polynomial is the minimal familiar two-root Mellin pole annihilator augmented by an independent linear phase factor; with three bumps the two constraints leave only one complex direction, whereas four allow a 2D pole-null family. The above identities have been **independently checked in exact symbolic algebra** in this investigation; the interval separation is an elementary analytic existence argument, not an already compiled Lean source.

### Exact prime computation and scope

Let \(c_0=1,c_1=\beta-s,c_2=1-s\beta,c_3=\beta\). The only active two-orientation cross-correlations occur at 3,9,27. Direct expansion gives, with \(x=\operatorname{Re}\beta\), \(\kappa_{3^k}=\log3/3^{k/2}\),
\[\frac{P_{\rm original}(v_\beta,v_\beta)}J=2\sum_{k=1}^3\kappa_{3^k}\operatorname{Re}\sum_{j=0}^{3-k}\overline{c_j}c_{j+k}=\log3\left(-\frac{28}{3}+4\sqrt3\,x\right).\]
Dividing by \(N(\beta)\) yields the stated ratio. The source support is inside \((-a,a)\) for \(a>\tfrac32\log3\), and the number of active von Mangoldt shifts is finite and exactly controlled by narrowing the bump. The full Gamma contribution generally differs among these tests, and may compensate a prime-term difference. Therefore **termwise** phase discrimination is proved; **full quadratic-value** discrimination is not claimed here.

**Important scope distinction:** the \(\beta=\pm1\) sources are real and satisfy the 2026 Theorem 8.1's real/pole-null input requirements. The \(\beta=i\) source is complex; extending that trace identity to complex tests requires a separate polarization/trace-class proof. Its scalar Weil explicit-formula prime identity remains valid. Thus there is a genuine **two-source real** trace-domain comparison already, not merely a formal imaginary-phase artifact.

## 6 — A true semilocal Hilbert map exists, but it distorts the native metric

Connes–Consani–Moscovici (2024), §4.6 eq. (57) and Theorem 4.6, constructs a **bounded invertible Hilbert-space isomorphism on Sonin spaces** \(\theta_S:\mathbf S_\lambda(\mathbb R,e_\infty)\to\mathbf S_\lambda(X_S,\alpha)\). In Mellin coordinates it is multiplication by the explicit **nonunit-modulus local factor**
\[D_S(t)=\prod_{p\in S\setminus\{\infty\}}(1-p^{-1/2-it}),\qquad\mathcal F_\mu w_S(\theta_S f)(t)=D_S(t)\,\mathcal F_\mu w_\infty f(t).\]
Consequently, on the relevant Hilbert domains, Parseval gives
\[\boxed{\prod_{p\in S_f}(1-p^{-1/2})\|f\|_2\le\|\theta_S f\|_2\le\prod_{p\in S_f}(1+p^{-1/2})\|f\|_2,}\]
where \(S_f=S\setminus\{\infty\}\). Both constants are explicit and positive **for each finite set**. For \(S=\{\infty,3\}\), they are \(1\mp3^{-1/2}\); the map is generally **not an isometry in the native \(L^2\) norm**. Nor are these bounds uniform as the place set grows. Theorem 4.6's word “Hilbertian isomorphism” must not be silently strengthened to 'unitary'.

**Why this matters to our previous Morishita result:** CR-MO-001 constructed a **unitary on each selected prime orbit with intrinsic normalized Haar measure**. The CCM \(\theta_S\) is a **bounded isomorphism of full Sonin spaces with an explicit nontrivial spectral multiplier**. These maps have different carriers, domains and metrics; there is no proved commuting-square identifying Morishita's dynamical map with \(\theta_S\). This gives a **measurable metric-distortion test**, not a new equivalence.

## 7 — Three kernels and two noninterchangeable enlargement parameters

Keep the following objects *typed*: (i) \(Q_{\rm Son}(g)=\operatorname{Tr}(\vartheta(g)S\vartheta(g)^*)=\|\vartheta(g)S\|_{HS}^2\ge0\), when finite, has **equality condition \(\vartheta(g)S=0\)** in the *test/operator* domain; (ii) Sonin's space is \(\operatorname{Ran}S\), and in the 2020 quasi-inner paper is identified with \(\ker((\rho_\infty)_{22})\) in a **different Hilbert factor/block**; (iii) an off-diagonal quasi-inner block has a different kernel/range type. These are **not the same nullspace** (retained September [PRED-SEAM-17](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_17_SONIN_EQUALITY_20260914.json)).

Likewise **adding finite places** \(S\subset S'\) changes \(D_S\), \(\theta_S\), and the Sonin realization, while **enlarging physical support** \(a\mapsto b\) changes the admissible test space and can activate new prime-power translations. These operations have no established equality or commutation theorem. If \(p>e^{2a}\), the original explicit formula on tests supported in \((-a,a)\) has **zero p-power contribution**, even though adding p to the place set alters the spectral multiplier \(D_S\). The sum of operator terms must be handled globally before inferring any form-preservation or support-naturality.

## 8 — The November 2026 conditional Schur route, and exactly what is missing

Andrews–Connes–La Gatta (Nov 9, 2026), Theorem 8.1, proves the genuine **semilocal trace identity** on real smooth pole-null window tests:
\[QW_\lambda(f,f)=\operatorname{Tr}\big(\theta(f)(1-P-\widehat P)\theta(f)^*\big).\]
This is **not** simply the positive Sonin trace: \(1-P-\widehat P\) has negative pair sectors, which is why positivity remains difficult. Their Theorem 8.14 gives a strictly positive lower margin on an explicit finite-codimensional **conditioned test subspace** only **assuming hypothesis (H2')** in the semilocal setting. (H2') controls the high-index negative-pair trace tail uniformly by a logarithmic Sobolev form; the paper proves an archimedean instance but not the general semilocal estimate. Under positive control of the complement, Lemma 8.16 gives the exact Schur reduction \(\mathfrak S=A-BD^{-1}B^*\).

**The conditional re-entry target is genuinely sharp:** Proposition 8.18 requires **(H2'), (R1), (R2)** plus a transition-and-two-pole Schur lower bound \(\lambda_{\min}(\mathfrak S|_{V'})\ge-\eta(\lambda)\), \(\eta(\lambda)\to0\), for all sufficiently large \(\lambda\). Only **under those assumptions** does monotonicity imply the global positivity statement equivalent to RH. No current Weil-Lab F40 gate count or selected two-bump Gram certificate supplies H2', R1/R2, or the uniform Schur estimate.

## 9 — Return verdict / next independent author

**PASS — exact source dictionary:** smooth original signed Weil form ↔ multiplicative-window formula via \(u=e^x\), with both poles, Gamma and full von Mangoldt correlations.

**PASS — falsifiable pole-null calibration:** two-bump IV030 span has **trivial two-pole-null intersection**, and the explicit four-cell Mellin-annihilator family repairs the admissibility while preserving same support/mass/prime interaction graph and a **three-valued genuine prime phase fingerprint**. Only the two real members are immediately in the published real trace-domain theorem.

**PASS — bounded metric transport:** CCM finite-S Sonin isomorphism has explicit product distortion bounds; **not a global Weil-form unitary**.

**OPEN — substantive new theorem:** global/larger support native spectral coercivity, semilocal H2' and R1/R2, uniform transition/pole Schur control, a true Morishita→CCM or CCM→original-RH selected-Picard carrier intertwiner. Original STRICT and BAL STOP records remain untouched. Nothing here establishes RH or a new general-form positivity theorem.

**Next: CR-DE-001 — Christopher Deninger.** Compare the conjectural arithmetic cohomological pairing and flow scaling with the actual metric distortion and phase/quantifier constraints above. Do not promote Morishita's flow reversal or CCM's Hilbertian Sonin correspondence into a Hodge star or a positive cohomological form without constructing it.
