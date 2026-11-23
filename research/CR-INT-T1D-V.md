# CR-INT-T1D-V — Two genuine broad-support Ritz directions per parity, with exact H¹ Gram certificates

**Canonical Return · November 23, 2026 · independent source-level analytical continuation of [T1D-U](CR-INT-T1D-U.md).**  
**Disposition:** **DERIVED / ANALYTIC** finite-multiplicity *lower* certificates for the **original** parity-separated compact Birman–Schwinger operator, supported by an executable exact-rational/SymPy audit. **Not** a rigorous evaluation of the full mollified Gamma or Schur **matrix entries**, not a certified leading-eigenvalue **upper** bound, not the spectral high/low coupling `eta_N`, not a proof of full complement positivity, semilocal domain equivalence, Lean, or RH. Draft PR #2 remains open and unmerged.

## 1. Freeze the original form, the metric, and the real objective

Use Suzuki's original signed Weil form on the unchanged `I=(-2,2)` and genuine compact smooth two-Mellin-pole-null source

\[
 \mathcal D^0=\{f\in C_c^\infty(I;\mathbf C):\ \int_I f(x)e^{\pm x/2}\,dx=0\}.
\]

The physical Hilbert space is `H=L²(I,dx)`. The original T1C unit selected directions `E=span{u_+,u_-}` and two pole carriers `P=span{e^{+x/2},e^{-x/2}}` remain fixed. Write `H_F=(E+P)^perp`, `F_s=D^0 ∩ H_F`. The complete source Weil prime term is `-<Af,f>`, where `A` is the original compressed 24-prime-power arithmetic operator of T1D-Q/U, and `gamma` is the full original Gamma Fourier-multiplier form. The genuine selected-source Schur correction `R` of T1D-E is positive semidefinite, with **no assumption that it vanishes** for the new tests.

The T1D-R operator is

\[
 D_F=\Gamma_F+6I\ge\tfrac58I,\qquad
 \mathscr B_F=D_F^{-1/2}(6I+A_F+\mathcal R)D_F^{-1/2},
\]

and its source-level equivalence is `q|F_s-R>=0` if and only if `sup sigma(B_F)<=1`. It is compact and self-adjoint, splitting into even and odd parity blocks. We seek *actual* source-domain finite-rank information—not the fictional inference that a positive compression bounds the operator from above. Suzuki's underlying finite-window compactness is prior canonical mathematics; the physical Schur-specific operator is our analytical construction, not a Suzuki theorem or a Connes Sonin intertwiner.

## 2. Four independent authentic source functions

Let `a=1/4`, `eps=1/10000`, and fix **any single even nonnegative C∞ probability mollifier** `rho_eps` with support in `[-eps,eps]`. For `sigma=e,o` put `p_e(x)=1`, `p_o(x)=x`, and for `k=3,4` define

\[
 b_{\sigma,k}(x)=p_\sigma(x)(1-16x^2)^k\mathbf 1_{|x|<1/4},\quad
 y_{\sigma,k}=(\partial_x^2-\tfrac14)(b_{\sigma,k}*\rho_\varepsilon).
\tag{V.1}
\]

Every `b` belongs to `H³` and is supported in `[-a,a]`; its derivatives to order two have zero boundary traces because `k>=3`. The `y` are nonzero C∞ compact tests supported in `[-a-eps,a+eps]`. Integration by parts (twice) gives **both** source moments exactly:

\[
 \int y_{\sigma,k}(x)e^{\lambda x}\,dx
 =(\lambda^2-\tfrac14)\int (b_{\sigma,k}*\rho_\varepsilon)(x)e^{\lambda x}\,dx=0
\quad (\lambda=\pm\tfrac12).
\tag{V.2}
\]

All four tests are physically orthogonal to the original E: the smallest absolute support center of an original selected bump is `log3/2 > 1/2`, separated from `a+eps=0.2501`. The even and odd subspaces are invariant under reflection. The two `k=3,4` functions in each parity are linearly independent: **the exact raw physical Gram is strictly positive definite**, and Young's inequality plus the uniform derivative-bound argument below show their common smoothing map is injective on each two-dimensional span.

The support diameter `1/2+2eps=0.5002 < 2/3 < log2` is **strictly smaller than every genuine prime-power translation**, yielding the operator identity

\[
 P_J A P_J=0,\qquad J=(-a-\varepsilon,a+\varepsilon).
\tag{V.3}
\]

Thus **every** pairing within the entire four-dimensional physical test space has **zero prime component** (not merely its diagonal): `q(y_i,y_j)=gamma(y_i,y_j)`. All pole contributions also vanish. This is a true complete-original-source simplification, not a sign change or truncation. However, the selected-to-central couplings `q(u_sigma,y_i)` in the Schur penalty need **not** vanish: their exact matrix remains uncomputed.

## 3. Exact two-by-two raw mass and derivative Gram matrices

Write `y^0_{sigma,k}=(d²/dx²-1/4)b_{sigma,k}` before mollification. Define the **rational** matrices

\[
 M_\sigma=(\langle y^0_{\sigma,i},y^0_{\sigma,j}\rangle)_{i,j=3,4},\quad
 K_\sigma=(\langle(y^0_{\sigma,i})',(y^0_{\sigma,j})'\rangle)_{i,j=3,4}.
\tag{V.4}
\]

Every entry is evaluated **exactly** by SymPy rational integration of polynomials on `[-1/4,1/4]`; no Fourier grid, floating quadrature, or unsmoothed test is substituted in the final domain. The checker obtains:

\[
 M_e=\begin{pmatrix}
 28194976/15015&92258752/45045\\
 92258752/45045&643216384/255255
 \end{pmatrix},\quad
 K_e=\begin{pmatrix}
 259883264/385&9012127744/15015\\
 9012127744/15015&3024257024/3465
 \end{pmatrix},
\]

\[
 M_o=\begin{pmatrix}
 184586/3465&5379076/109395\\
 5379076/109395&762986944/14549535
 \end{pmatrix},\quad
 K_o=\begin{pmatrix}
 23448144/715&358195328/15015\\
 358195328/15015&1498703872/51051
 \end{pmatrix}.
\tag{V.5}
\]

The checker verifies `M_sigma>0` and **strict rational Sylvester criteria**

\[
 \boxed{37^2M_e-K_e>0,\qquad45^2M_o-K_o>0}.\tag{V.6}
\]

These inequalities hold for **every** real or complex coefficient pair, by complexification of the real symmetric Gram. For any `y0` in the corresponding raw span, write `y=y0*rho_eps`. Positivity and normalization of `rho_eps` imply `||y'||<=||y0'||`; Minkowski and the fundamental theorem of calculus imply `||y-y0||<=eps||y0'||`. Hence

\[
 \frac{\|y'_e\|}{\|y_e\|}<\frac{37}{1-37/10000}<38,\qquad
 \frac{\|y'_o\|}{\|y_o\|}<\frac{45}{1-45/10000}<46
\tag{V.7}
\]

for **every nonzero vector of the two-dimensional smooth** even and odd subspaces. In particular the smoothing map does not collapse a nonzero raw combination.

## 4. Every direction in the two-dimensional spaces gives a genuine compact Ritz lower witness

Import T1D-U's separately audited **analytical** digamma/Jensen inequality, including `m_Gamma(0)<-5372/1000`. For `V=||y'||²/||y||²`, its explicit rational upper formula at `M=50` is

\[
 \gamma[y]/\|y\|²< -\tfrac{5372}{1000}
 +\sum_{n=0}^{50}\frac{V}{(n+1/4)(4(n+1/4)^2+V)}
 +\tfrac12\log\bigg(1+\frac{V}{4(50+1/4)^2}\bigg).
\tag{V.8}
\]

The sum and logarithmic tail bound are increasing in `V`, so (V.7), the exact-fraction finite sum and an atanh series rational *upper* log enclosure give uniformly on the entire respective smooth two-spaces

\[
 \boxed{\gamma[y_e]/\|y_e\|²<181/100,\qquad
 \gamma[y_o]/\|y_o\|²<2.}
\tag{V.9}
\]

The checker prints the rationally bounded Jensen-upper numerical diagnostics `~1.801096919749` and `~1.992636390503`. These are values of an **upper expression**, not actual Gamma energies.

For `0 != y` in either two-space, put `z=D_F^{1/2} y`. Since the complete original source prime self-pairing on `J` is zero, and since the **genuine Schur term is positive**, the exact Rayleigh quotient satisfies

\[
 \frac{\langle \mathscr B_F z,z\rangle}{\|z\|²}
 =\frac{6\|y\|²+\langle\mathcal R y,y\rangle}{6\|y\|²+\gamma[y]}
 \ge\frac6{6+\gamma[y]/\|y\|²}.
\tag{V.10}
\]

By (V.9), the **Courant–Fischer min–max principle** now gives **two genuine positive eigenvalues in EACH parity sector** with lower bounds

\[
 \boxed{\lambda_{2,+}(\mathscr B_F)>\frac{600}{781}=0.768245\ldots,\qquad
 \lambda_{2,-}(\mathscr B_F)>\frac34=0.75.}
\tag{V.11}
\]

Here `lambda_2` denotes the second-largest positive eigenvalue of the compact self-adjoint parity block; (V.11) counts multiplicity. **These are lower, not upper, spectral bounds.** They do not give an eigenbasis or any guarantee about the entire remaining complement.

## 5. Better *leading* lower witnesses by exact rational polynomial mixing

The two raw profiles can be mixed using fixed rational coefficients chosen from an exploratory generalized derivative Rayleigh calculation, but the **actual certificates are exact rational checks**, not floating-point optimization:

\[
 y_e^*=y_{e,3}+\tfrac{13}{12}y_{e,4},\qquad
 y_o^*=y_{o,3}+\tfrac95 y_{o,4}.
\tag{V.12}
\]

The raw exact Gram matrices (V.5) certify

\[
 \frac{\|(y_e^0)^*{}'\|}{\|(y_e^0)^*\|}<18,
 \qquad
 \frac{\|(y_o^0)^*{}'\|}{\|(y_o^0)^*\|}<\frac{578}{25}=23.12.
\]

After the **same** fixed `eps=1/10000` mollification, they give `||y_e^{*'}||/||y_e^*||<181/10` and `||y_o^{*'}||/||y_o^*||<116/5`. The same rigorous U-type Gamma bound (V.8) verifies

\[
 \boxed{\gamma[y_e^*]/\|y_e^*\|²<53/50,
 \qquad\gamma[y_o^*]/\|y_o^*\|²<131/100.}
\tag{V.13}
\]

Hence the **actual** parity compact spectral maxima satisfy

\[
 \boxed{\sup\sigma(\mathscr B_{F,+})>\frac{300}{353}\approx0.84985835,\qquad
 \sup\sigma(\mathscr B_{F,-})>\frac{600}{731}\approx0.82079343.}
\tag{V.14}
\]

This improves T1D-U's `>60/73≈0.821918` (even) and `>300/377≈0.795756` (odd), and is independent of T1D-I's controversial numerical two-sided residual bounds: only the undisputed **positive Schur form sign** is used. The exact gamma and R cross-Gram entries remain **uncomputed**. No claim that the leading eigenvalues are close to, or below, 1 follows.

## 6. Full formal finite matrix, missing data, and stopping obligations

For `i,j=3,4` in parity `sigma`, the genuine full-source **form** numerator and denominator matrices on these smooth tests are

\[
 C_{ij}=6\langle y_i,y_j\rangle
 +\sum_{\tau=\pm}\frac{q(y_i,u_\tau)\overline{q(y_j,u_\tau)}}{\mu_\tau},
\quad
 D_{ij}=6\langle y_i,y_j\rangle+\gamma(y_i,y_j).
\tag{V.15}
\]

The prime quadratic matrix **vanishes exactly**, because `P_J A P_J=0`, while the finite rank Schur terms are **retained**. The Rayleigh quotient (V.10) is the generalized pencil `C v=lambda D v` on the two-dimensional source test space. (V.5) gives **exact pre-mollification mass and derivative Gram matrices**, *not* the complete mollified `C,D` matrices; (V.9) and (V.13) certify one-sided bounds without inventing their unknown entries. Exact full mollified Gamma integrals and all source cross-terms `q(u_tau,y_i)` still must be computed or independently enclosed.

**Source/audit custody:** `src/audits/canonical_return_t1d_v_central_ritz_multiplicity.py` checks the exact rational Gram integrals, positive definiteness, sharp rational mixed-vector inequalities, mollification robustness, and finite-sum Jensen budgets. It imports T1D-U's analytic digamma bound and its separate certified `m_Gamma(0)` rational enclosure, as well as the established T1D-R compact source form. It does **not** calculate actual eigenvalues of `D_F`, `eta_N`, the actual full two-by-two Gamma/Schur matrices, or any top spectral **upper** enclosure. No Lean kernel or independent peer review.

**Next CR-INT-T1D-W:** attempt rigorous **two-by-two source matrix intervals** for the genuinely mollified profiles (especially the off-diagonal Gamma integral), and real selected–central Schur cross terms. In parallel, seek a source-valid control on `eta_N` or an alternate high-energy norm, with no inference from finite Ritz *lower* witnesses to the unresolved `sup sigma(B_F)<=1`. The existence of *four* well-evidenced positive compact spectral modes does **not** establish the sign of the original full Weil complement or any Phase IV theorem; the full RH claim remains **OPEN / NO TRANSFER**.