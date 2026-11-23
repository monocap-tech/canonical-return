# CR-INT-T1D-T — A structure-adapted quadratic Schur tail and broad-support pole-null Ritz probes

**Canonical Return · 2026-11-23 (America/Los_Angeles) · independent smooth-source analysis.**  
**Disposition:** **ANALYTICAL / CONDITIONAL CERTIFICATION RULE** with exact-rational/SymPy audit of finite constants and concrete source probes. **No top Ritz eigenvalue has been computed; no finite low-mode Schur matrix is certified; no whole-complement positivity, completed Suzuki–Connes correspondence, Lean theorem or RH proof.** This is the continuation of [T1D-R](CR-INT-T1D-R.md) and [T1D-S](CR-INT-T1D-S.md), in the existing draft PR #2.

## 1. Preserve the exact source, domain and earlier ownership

Use Suzuki's **original complete signed Weil form** on compact smooth logarithmic tests in `I=(-2,2)`, including the two zeta pole terms, negative von Mangoldt prime-power sum and Gamma distribution ([arXiv:2606.09096v4 §1.1 and §2.4](https://arxiv.org/html/2606.09096v4)). Fix the original width `ell=1/10000` T1C source and its physical-orthonormal pole-null pair `u_+`, `u_-`, whose selected Weil values are `mu_±>0`. The smooth pole-null core and its complementary physical Hilbert space are

\[
\mathcal D^0=\{f\in C_c^\infty(I):\int_I f(x)e^{\pm x/2}dx=0\},\quad
\mathscr H_F=(\operatorname{span}\{u_+,u_-,e^{x/2},e^{-x/2}\})^\perp,
\quad F_s=\mathcal D^0\cap\mathscr H_F.
\]

T1D-E–I define the exact finite-rank **positive** Weil cross-block penalty `R` on this physical complement. The actual signed quadratic Schur remainder is `q(y,y)-<Ry,y>`. Set `A` to the original **positive-kernel** (but not positive-semidefinite) joint prime translation operator; the arithmetic part of `q` is `-<Ay,y>`. By T1D-R's specialization of Suzuki's **canonical** finite-window compactness theorem, the strictly positive shifted Gamma form `d[y]=gamma[y]+6||y||²` defines `D_F>=5/8` with compact resolvent on `H_F`. The exact compact Birman–Schwinger operator is

\[
\boxed{\mathscr B_F=D_F^{-1/2}(6I+A_F+\mathcal R)D_F^{-1/2},\quad A_F=\Pi_F A\Pi_F.}\tag{T.1}
\]

The unresolved, unchanged complete smooth-core Schur condition is **exactly** `sup sigma(B_F)<=1`. The reflection parity sectors are invariant. **No equivalence with the Connes/Sonin conditioned operator is assumed**; see [Andrews–Connes–La Gatta, arXiv:2611.13129v1 §8.6](https://arxiv.org/html/2611.13129v1).

## 2. The exact improvement: the `6I` term cannot mix Gamma eigenmodes

Let `P_N` be the spectral projection onto the first `N` eigenvectors of `D_F`, with `Q_N=I-P_N`. Write `h=d_{N+1}>0` for the first omitted Gamma eigenvalue, and let `d_1>=5/8` be the least Gamma eigenvalue. Decompose `B_F` into the `P_N` and `Q_N` blocks:

\[
\mathscr B_F=\begin{pmatrix}B_N&X_N^*\\ X_N&T_N\end{pmatrix}.
\]

Crucially, `P_N` is a **spectral** projection of `D_F`. The term `6D_F^{-1}` commutes with it and makes **no** off-diagonal contribution:

\[
\boxed{X_N=Q_ND_F^{-1/2}(A_F+\mathcal R)D_F^{-1/2}P_N.}\tag{T.2}
\]

Write `L=||A_F+R||`, with `A_F` retaining its original signed quadratic role in (T.1). Then

\[
\boxed{\|T_N\|\le\frac{6+L}{h},\qquad
\|X_N\|\le\frac{L}{\sqrt{d_1h}}.}\tag{T.3}
\]

These are **source-valid operator bounds**; (T.2) is stronger than naively using `||6I+A_F+R||` in the off-diagonal block. No claim of extra smoothing of `A` is made.

### Exact block Schur complement

If `h>6+L`, then `I-T_N` is strictly positive. By completing the square **in the two spectral blocks**,

\[
\boxed{I-\mathscr B_F\ge0\quad\Longleftrightarrow\quad
I-B_N-X_N^*(I-T_N)^{-1}X_N\ge0.}\tag{T.4}
\]

A rigorous **sufficient**, not necessary, scalar stopping rule is therefore

\[
\boxed{\lambda_{\max}(B_N)+\frac{L^2}{d_1(h-6-L)}\le1
\quad\Longrightarrow\quad q(y,y)\ge\langle\mathcal Ry,y\rangle\ \forall y\in F_s.}\tag{T.5}
\]

The new tail penalty decreases as `O(1/h)` rather than T1D-R's conservative *additive operator-norm approximation* `O(h^(-1/2))`. This is a **quadratic Schur correction**, not a claim that `||B_F-P_NB_FP_N||` itself decays as `O(1/h)`. It is essential not to confuse those two statements.

## 3. The new finite-dimensional diagnostic: an exact residual Gram

For a concrete `D_F`-orthonormal eigenbasis `e_1,...,e_N` on the physical complement, let

\[
\boxed{\eta_N^2:=\|Q_N(A_F+\mathcal R)D_F^{-1/2}P_N\|^2.}\tag{T.6}
\]

It is the top eigenvalue of an **N×N positive Gram matrix** built from the genuine arithmetic-plus-Schur images:

\[
\boxed{(G_N)_{ij}=\frac{\langle Q_N(A_F+\mathcal R)e_i,\ Q_N(A_F+\mathcal R)e_j\rangle}{\sqrt{d_id_j}}.}\tag{T.7}
\]

Then `||X_N||²<=eta_N²/h`, giving the **sharper form-adapted sufficient gate**

\[
\boxed{\lambda_{\max}(B_N)+\frac{\eta_N^2}{h-(6+L)}\le1.}\tag{T.8}
\]

All entries of (T.7) are **original-source physical pairings**; one can obtain them from the `N` images of `(A_F+R)e_j` and their projections, without choosing an infinite complementary basis. However, this **requires actual** `D_F` eigenvectors or a rigorously equivalent projection; neither `eta_N` nor `lambda_max(B_N)` has yet been evaluated. Equation (T.8) is a falsifiable certification target, not a positivity result.

### Conservative numerical improvement and its limitation

Import the separately audited **T1D-Q** estimate `||A||<203/25`, and **conditionally** import T1D-I's still-unreviewed `||R||<155/182`. Then

\[
L\le\|A_F\|+\|\mathcal R\|<9,\quad 6+L<15,\quad d_1\ge5/8.
\]

For `h>15` the fully explicit sufficient condition becomes

\[
\boxed{\lambda_{\max}(B_N)+\frac{648}{5(h-15)}\le1.}\tag{T.9}
\]

At `h>=275` the new tail allowance is `<162/325<1/2`. Compare T1D-S: the prior *specific* conservative Ritz-norm estimator could not be below `1/2` without `h>=5760`. **This is a reduction in the required Gamma-energy threshold of a sufficient estimator, not a certified reduction to 275 computable modes.** Since `d_N` grows only logarithmically, even attaining `h=275` may require an astronomically large rank. From T1D-S's min–max bound `d_{N+1}<=46+log(1+pi(N+5)/4)`, having `h>=275` is possible only if `N >= (4/pi)(exp(229)-1)-5` (with the inequality interpreted as an exact necessary lower bound). No practical low-mode certificate follows.

Without T1D-I, (T.4)–(T.8) remain **unconditional algebraic/form criteria**, but the numeric constants `9`, `15`, `648/5` and the comparison at 275 are **conditional**. This preserves the established evidence-transfer boundary.

## 4. An explicit pair of broad-support, source-admissible parity probes

To avoid examining only the highly localized T1D-J/K witnesses, construct two genuine smooth **moderate-frequency** tests *independently* of the selected carrier's `ell=10^-4`:

\[
\begin{gathered}
b_e(x)=(1-16x^2)^5\mathbf 1_{|x|<1/4},\qquad
b_o(x)=x(1-16x^2)^5\mathbf 1_{|x|<1/4},\\
\rho\in C_c^\infty(-1/200,1/200),\ \rho\ge0,\ \rho(-x)=\rho(x),\ \int\rho=1,\\
\boxed{y_e=(\partial_x^2-1/4)(b_e*\rho),\quad
       y_o=(\partial_x^2-1/4)(b_o*\rho).}
\end{gathered}\tag{T.10}
\]

They are nonzero, smooth, compactly supported in `[-51/200,51/200]`, and even/odd respectively. For `lambda=+1/2` and `lambda=-1/2`, integration by parts on the compact smooth source gives

\[
\int y_\sigma(x)e^{\lambda x}dx=(\lambda^2-1/4)\int(b_\sigma*\rho)(x)e^{\lambda x}dx=0.
\]

Thus each test is **exactly two-Mellin-pole-null**. The nearest four-cell selected center has absolute coordinate `log3/2>1/2`; its support width is at most `10^-4`, so the support of each new `y` is physically disjoint from `E`. Therefore `y_e,y_o\in F_s`, with opposite parity and physical orthogonality.

**Entire original prime contribution vanishes for each of these self-pairings:** their support diameter is `51/100<2/3<log2`, whereas every genuine von Mangoldt shift has magnitude at least `log2`. The two zeta pole terms vanish as above. Hence

\[
\boxed{q(y_\sigma,y_\sigma)=\gamma[y_\sigma]\quad(\sigma=e,o),}\tag{T.11}
\]

with **no** prime omission or substitution of another form. These probes do **not** by themselves show Gamma positivity or provide a negative Weil direction.

The companion exact SymPy checker integrates the unmollified polynomial energy exactly:

| Probe | `||y_0||²` | `||y_0'||²` |
|---|---|---|
| even | `9581387776/2909907` | `3906493153280/2909907` |
| odd | `3621357056/66927861` | `100490301440/2909907` |

Here `y_0=(partial²-1/4)b_σ`; its `H¹` norm is finite and it is nonzero. Young/Minkowski yield `||y-y_0||₂<=||y_0'||₂/200` and `||y'||₂<=||y_0'||₂`, so the exact rational inequalities imply

\[
\boxed{\|y'_e\|_2/\|y_e\|_2<24,\qquad
\|y'_o\|_2/\|y_o\|_2<30.}\tag{T.12}
\]

These are **broad-support compared to the selected 10^-4 carrier**, with explicit finite Fourier-frequency first-moment control; they are not bandlimited. As a deliberately loose diagnostic, T1D-F's global Gamma upper multiplier `m_Gamma(xi)<=log(1+|xi|)+40`, Plancherel and Jensen imply `gamma[y]/||y||²<40+log31<44`. Therefore `d[y]/||y||²<50`. Since `6+<Ay,y>/||y||²+<Ry,y>/||y||²>=6`, the actual compact parity Birman–Schwinger spectra have **independent** trial lower bounds

\[
\boxed{\sup\sigma(\mathscr B_{F,+})>3/25,\qquad
\sup\sigma(\mathscr B_{F,-})>3/25.}\tag{T.13}
\]

These are weaker than T1D-S's `>7/25` and `>31/100`, but probe **qualitatively different, moderate-frequency** directions. No upper bound below the critical 1 is inferred. The exact self Gamma values for these tests remain to be evaluated if stronger spectral information is desired.

## 5. Verification, falsifiers and the next gate

**Established analytically:** exact low/high spectral Schur identity (T.4), no off-diagonal coupling from `6I`, quadratic error certificate (T.5), finite Gram for mixed operator (T.7) and source-preserving moderate-frequency parity probes (T.10)–(T.12). **Exact arithmetic/symbolic audit PASS:** the conditional constant derivation at `h=275`, finite 2×2 Schur algebra and exact broad polynomial derivative energies. The audit does **not** derive Suzuki's Gamma multiplier, the compact-resolvent theorem, the complete Weil representation, nor the actual spectrum or any numerical `eta_N`.

**Falsifiers/remaining tasks:** compute/interval-enclose `eta_N`, `lambda_max(B_N)` and a genuine **lower** bound for `d_{N+1}` in the same spectral domain; independent review of imported T1D-I norm estimates; compute true Gamma values for the broad `y_e,y_o`; determine whether the low–high mixed prime block has additional decay; add the original pole directions and any support-independent estimates; prove a proper Suzuki–Connes semilocal form/metric/domain dictionary. **No RH claim.**

**Next CR-INT-T1D-U:** investigate the arithmetic off-diagonal *directly*: exact compressed-shift images of a computable low-mode test basis, plus an independently certified lower bound on the complementary shifted Gamma energy. Without those two quantities, the improved Schur criterion does not close. In parallel, evaluate the actual source Gamma integrals of `y_e,y_o` to decide whether moderate-frequency trials move the leading Ritz values nearer the threshold 1. No transfer to Phase IV until full source/metric/domain preservation and a relevant theorem are verified.

**Reproduction:** `python3 src/audits/canonical_return_t1d_t_tail_schur.py` (requires SymPy). This is **not** a Lean or CI certificate.