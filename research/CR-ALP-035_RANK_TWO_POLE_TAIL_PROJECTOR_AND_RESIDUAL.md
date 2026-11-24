# CR-ALP-035 — Exact pole-null high-sine projector and the original signed Weil Schur residual

**2026-11-23 PST; Canonical-Return analytical research; DRAFT/NONPROMOTING.** This is a new exact geometric reduction of the missing infinite-dimensional Schur coupling, plus exact-rational arithmetic/finite-algebra audit. The true 15,000-dimensional original signed finite Schur inequality, its spectral acceptance, the full fixed-collar positivity and RH remain OPEN. None of the continuum arguments in this note has been newly Lean-compiled.

## 0. Source custody and retained results

The original Zeta23 source is still `monocap-tech/weil-lab` branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`; its `literatureRHS` consists of **two original poles, minus paired von Mangoldt terms, plus the full Gamma/digamma bracket**. The prime support at the fixed physical radius `a=11/10` consists of the same seven powers `{2,3,4,5,7,8,9}`. The original two-pole-null function space is

\[
\mathscr H_0=\{f\in L^2(-a,a;\mathbb C):\langle f,\cosh(x/2)\rangle=\langle f,\sinh(x/2)\rangle=0\}.
\]

As in CR-ALP-025, define the orthonormal Dirichlet sine basis on `I=(-a,a)` by `e_j(x)=a^-1/2 sin(j π(x+a)/(2a))`. Let `S_N` project onto `e_1,...,e_N`, `R_N=I−S_N`, and `Π_0` project onto `H_0`. The physical finite space is `E_N=Π_0 ran(S_N)`; its physical orthogonal complement in `H_0` is denoted `Q_N H_0`.

CR-ALP-027 proves analytically, not in Lean, that for `N=15,000`, all logarithmic form-domain `q` in `Q_N H_0` satisfy **the original signed Weil** tail inequality `W(q,q)≥(1/4)||q||²`. CR-ALP-034 certifies the independent twelve-mode lower bound `3×10^-12` in the true physical metric, but it does **not** certify `E_15000` or its cross-coupling to `Q_15000`.

## 1. NEW EXACT INFINITE-DIMENSIONAL PROJECTOR IDENTITY

Put `b_e(x)=cosh(x/2)`, `b_o(x)=sinh(x/2)`, and `B: C²→L²(I)` by `B(z_e,z_o)=z_e b_e+z_o b_o`. The two pole constraints are `B* q=0`.

**Lemma.** For every finite `N`,

\[
\boxed{Q_N\mathscr H_0=\operatorname{ran}(R_N)\cap\ker B^*.}\tag{35.1}
\]

Indeed `q∈H_0` is orthogonal to `Π_0e_j` precisely when it is orthogonal to `e_j`. Hence the usual sine-orthogonality constraints and the pole moments are independent conditions.

Define the 2×2 **high-frequency pole Gram**

\[
G_N=B^*R_NB.
\]

`G_N` is positive definite: the tail sine coefficients of `cosh(x/2)` are nonzero on infinitely many odd indices and those of `sinh(x/2)` are nonzero on infinitely many even indices. Reflection makes `G_N` diagonal. Consequently the actual orthogonal projector onto the whole (infinite-dimensional) tail is **exactly**

\[
\boxed{Q_N=R_N-R_NB\,G_N^{-1}B^*R_N.}\tag{35.2}
\]

No infinite Gram matrix is inverted; the subtraction has **rank exactly two** on the high-sine subspace. Formula (35.2) acts on any `h∈L²(I)`, not merely on an already pole-null vector.

## 2. NEW EXACT ORIGINAL-SIGNED CROSS-TAIL RESIDUAL FORMULA

Let `T` be the full original physical signed Gamma-minus-seven-prime operator on the fixed collar and let `A=Π_0 T Π_0` be its pole-null compression. The genuine poles vanish in mixed forms involving a pole-null source. For any `p∈E_N∩Dom(A)` write `h=A p`. Then

\[
\boxed{\|Q_NAp\|_2^2=\|R_Nh\|_2^2-(B^*R_Nh)^*G_N^{-1}(B^*R_Nh).}\tag{35.3}
\]

The same formula holds with `h=T p`, since `Q_NΠ_0=Q_N`, and using `T` avoids unnecessary pole projection inside high-sine coefficient calculations.

If `φ_i` are physical finite-section columns and `a_{ji}=⟨e_j,T φ_i⟩`, while `c_{αj}=⟨e_j,b_α⟩`, then the *exact positive semidefinite* infinite-tail residual matrix is

\[
\boxed{\mathcal B_{N,ik}=\sum_{j>N}\overline{a_{ji}}a_{jk}-d_i^*G_N^{-1}d_k,\quad
 d_i=\sum_{j>N}(\overline{c_{ej}}a_{ji},\overline{c_{oj}}a_{ji})^\mathsf T.}\tag{35.4}
\]

The sums converge by Bessel/Parseval because the fixed physical source columns lie in `Dom(A)` (Fourier decay O(1/r) for the two boundary anchors suffices against the logarithmic multiplier). The exact Gram correction must be evaluated with the **actual physical metric**; it is not an identity-matrix surrogate.

By parity, for `p` of pure reflection parity the corresponding single `b` gives the **weighted-variance identity**

\[
\boxed{\|Q_NAp\|_2^2=\min_{t∈\mathbb C}\sum_{j>N,\ \mathrm{matching\ parity}}|a_j-tc_j|^2.}\tag{35.5}
\]

This cancellation is part of the actual original pole-null geometry. Dropping it gives a valid but potentially very wasteful upper bound. Neither (35.3) nor (35.5) asserts that the residual is small enough for Schur positivity.

## 3. A QUANTITATIVE LOW-TO-TAIL OPERATOR BOUND FOR ALL INTERIOR SOURCES

For any `h∈H¹(I)` choose the **unique** coefficients

\[
 t_e=\frac{h(a)+h(-a)}{2\cosh(a/2)},\qquad
 t_o=\frac{h(a)-h(-a)}{2\sinh(a/2)}
\]

and set `h_0=h−t_e b_e−t_o b_o`. By construction `h_0∈H¹_0(I)`. Since the true projector kills both pole functions, `Q_Nh=Q_N h_0`. The Dirichlet sine/derivative Parseval identity gives

\[
\boxed{\|Q_Nh\|_2 \le \|R_Nh_0\|_2
\le\frac{2a}{\pi(N+1)}\|h'−t_e b_e'−t_o b_o'\|_2.}\tag{35.6}
\]

It follows for **every** interior combination `p∈span{w_3,...,w_M}`, `M≤N`, that

\[
\boxed{\|Q_N A p\|_2^2\le \frac{4a^2}{\pi^2(N+1)^2}
\| (Tp)'−t_e(p)b_e'−t_o(p)b_o'\|_2^2.}\tag{35.7}
\]

**Domain justification.** Each original `w_j` is a finite Dirichlet sine combination with zero endpoint values. Its zero-extended Fourier transform is `O((1+|r|)^−2)` by two integrations by parts including endpoint derivative jumps. Since the *full original signed symbol* obeys `|m_a(r)|≤log(1+|r|)+16` (CR-ALP-028), `m_a\hat w_j` has an `L²` first Fourier derivative; hence `T w_j∈H¹(ℝ)` and its restriction lies in `H¹(I)`. Prime shifts preserve this regularity. No such `H¹` assertion is made for either exceptional boundary anchor `u_1,u_2`: their zero-extended Fourier tails are `O(1/r)` and the multiplier need not map them to `H¹`.

Thus the original interior-column cross-tail residual has a **rigorous computable derivative/trace upper form**. Let `φ_i` range over interior directions only and let `r_i=(Tφ_i)'−t_e(φ_i)b_e'−t_o(φ_i)b_o'`. Then

\[
\boxed{0\preceq\mathcal B_N^{\rm interior}\preceq
 \frac{4a^2}{\pi^2(N+1)^2}(\langle r_i,r_j\rangle)_{ij}.}\tag{35.8}
\]

This is a genuine analytical operator-norm control on `E_M^{interior}↔Q_N`. Its **numerical right-hand side has not been rigorously enclosed** for the original signed matrix; in particular no full 15,000-dimensional Schur criterion has passed.

The two anchor columns can be handled exactly by (35.3) or (35.5), or bounded via separate fractional Sobolev estimates; they **cannot** be inserted into the `H¹` derivative bound without a proof of the missing regularity.

## 4. SOURCE-INDEXED EXACT POLE GRAM AND RATIONAL ENCLOSURES AT N=15000

For odd `j`,

\[
\langle e_j,b_e\rangle=\frac{2\nu_j\cosh(a/2)}{\sqrt a\,(\nu_j²+1/4)},
\]

and all even-j `b_e` coefficients vanish. For even `j`,

\[
\langle e_j,b_o\rangle=-\frac{2\nu_j\sinh(a/2)}{\sqrt a\,(\nu_j²+1/4)},
\]

and all odd-j `b_o` coefficients vanish. Therefore

\[
\boxed{G_N=\begin{pmatrix}
\frac{4\cosh²(a/2)}a\sum_{j>N,\,j\ odd}\frac{\nu_j²}{(\nu_j²+1/4)²}&0\\
0&\frac{4\sinh²(a/2)}a\sum_{j>N,\,j\ even}\frac{\nu_j²}{(\nu_j²+1/4)²}
\end{pmatrix}.}\tag{35.9}
\]

An independent stdlib `Fraction` audit, using `3<π<22/7`, exact Taylor enclosures for hyperbolics at `11/20`, `ν_{N+1}>1`, the elementary comparison `(16/25)/ν²≤ν²/(ν²+1/4)²≤1/ν²`, and the parity-step integral test, gives strict rational envelopes for `N=15000`:

\[
\boxed{\frac1{20000}< (G_N)_{ee}<\frac1{10000},\qquad
\frac1{90000}<(G_N)_{oo}<\frac1{40000}.}\tag{35.10}
\]

The sharper exact rational endpoints (and decimal displays) are computed in the executable audit. In particular `G_N` is *quantitatively* invertible and `G_N^{-1}` is a small **explicit two-by-two** diagonal operator. These rough brackets are not a replacement for actual signed operator entries.

## 5. EXACT HIERARCHICAL SHELL PROJECTOR

The finite sections are nested, `E_M⊂E_N` for `M<N`, and their tail projectors satisfy `Q_N⪯Q_M` (as orthogonal projections). Consequently

\[
\boxed{\mathsf P_{M,N}=Q_M-Q_N}\tag{35.11}
\]

is the **exact physical orthogonal projector** onto `E_N⊖E_M`, hence positive semidefinite and of rank `N−M`. Formula (35.2) represents each `Q` using only ordinary sine tail sums and two-by-two pole Grams, enabling a hierarchical representation of every finite shell without a new large physical Gram inversion.

For example, the cross-tail energy of any `h∈L²(I)` splits exactly:

\[
\boxed{\|Q_Mh\|_2²=\|\mathsf P_{M,N}h\|_2²+\|Q_Nh\|_2².}\tag{35.12}
\]

This is a **residual decomposition**, not a claim that the corresponding signed Weil form on the intervening physical shell is positive. To certify full collar positivity one still needs the genuine original signed finite shell matrices and their Schur interactions.

## 6. RETAINED EXACT ACCEPTANCE TEST AND OPEN GATES

CR-ALP-027's analytical `W(q,q)≥1/4||q||²` on `Q_15000` yields the unchanged exact sufficient condition

\[
\boxed{H_{15000}-4\mathcal B_{15000}\succeq0\quad\Longrightarrow\quad
W_{\mathrm{original}}\ge0\ \text{on the complete two-pole-null fixed collar}.}\tag{35.13}
\]

Here `H` is the **genuine signed** finite physical form, and the newly derived `\mathcal B` is the **exact** rank-two-corrected infinite high-sine residual (35.4). Equations (35.2)–(35.12) supply a constructive, form-preserving basis for evaluating and bounding that residual, including finite shells. **They do not establish (35.13).** A certified positive 12-mode block from CR-ALP-034 neither decides the remaining `E_15000` directions nor the residual budget.

## 7. EXECUTABLE AUDIT AND GOVERNANCE

The standalone `src/audits/cr_alp_035_rank_two_high_sine_audit.py` is stdlib-only. **Locally executed PASS**: exact rational enclosures (35.10) at original radius `11/10` and finite cutoff `N=15000`; exact `Fraction` regression of idempotence, self-adjointness, pole and low-sine annihilation of the rank-two projector; exact Pythagorean residual downdate; and exact positivity/idempotence of a nested shell projector. Its finite-dimensional projection regression is an abstract toy, not a claimed finite zeta/Weil Gram. The analytic projector identity, its operator-domain hypotheses and derivative bound are **not machine-formalized in Lean**.

**CR-ALP-035 ledger:** SOURCE-PINNED original signed form; NEW infinite-dimensional **rank-two high-sine projection identity, analytical proof**; NEW exact residual as rank-two-corrected high-sine norm, analytical proof; NEW `H¹` interior cross-tail norm estimate, analytical proof; high-frequency pole Gram arithmetic PASS; true 15,000-mode matrix/residual computation OPEN; full collar PSD, support extension, F41 242 sufficiency, Suzuki–Connes identification and RH OPEN.

**Next CR-ALP-036:** compute **certified** sine coefficients and residual Gram entries for the actual original twelve source columns (especially the two exceptional anchors), using interval enclosures for the finite physical translation/digamma operator and analytic high-index envelopes. Evaluate the new exact `\mathcal B_{15000}` on the already-certified near-null twelve-mode directions and compare it to the physical `3×10^-12` margin. In parallel the 12→15000 signed shell must be certified. Do not promote a small-mode residual test alone to whole-collar positivity.

**Governance:** Keep Canonical-Return PR #3 DRAFT/NONPROMOTING. No Weil-Lab receiving changes, new Lean theorem, F41/F42 acceptance change, Hodge/provenance change, arbitrary-radius extension, or RH announcement.