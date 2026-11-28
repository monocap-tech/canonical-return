# CR-INT-T1D-K — Prime-4 translated probes, a six-dimensional Weil-positive source, and the residual left outside finite arithmetic probes

**Canonical Return · independent original-signed-Weil smooth-core research · November 23, 2026 (America/Los_Angeles).**  
**Status:** **ANALYTICAL RESTRICTED SOURCE RESULT; exact rational/SymPy audit PASS.** The finite-prime arithmetic, pole-null coefficients, physical Gram, rigorous mixed-pairing bounds and finite-dimensional Gershgorin margins are directly checked. The analytic Gamma source-integral bound and T1C's single-bump diagonal are separately stated inputs, not proved by SymPy. The quantitative *unmeasured-total-residual* lower bound additionally depends on T1D-I's **prior**, analytically derived two-sided finite-width residual bounds, which have not been Lean-accepted or independently peer-reviewed. No whole-complement theorem, Sonin/prolate subspace transport, global Weil positivity or RH claim.

**Source contract:** [Suzuki, *Weil’s quadratic form via the screw function*, arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4), original distribution including both poles, the complete negative von Mangoldt prime-power sum and the full Gamma distribution. Companion source notes: [T1C](../archive/space/10-t1c.md), [T1D-E](CR-INT-T1D-E.md), [T1D-H](CR-INT-T1D-H.md), [T1D-I](CR-INT-T1D-I.md), [T1D-J](CR-INT-T1D-J.md). The source-defined mixed representer is T1D-E's \(\mathcal K_2\), **not** Suzuki's inverse Neumann Laplacian. Connes semilocal domain and conditioning assumptions are **not** transported here.

## 1. Original source and the two new arithmetic probes

Fix \(I=(-2,2)\), \(\ell=1/10000\), \(d=\log3\), \(s=4/\sqrt3\), and the original even, nonnegative smooth compact bump \(\phi\) from T1C supported in an interval of width \(\ell\). Write \(J=\|\phi\|_2^2\), \(x_j=(j-3/2)d\), and \(\phi_j(x)=\phi(x-x_j)\). The old pole-null generators are
\[
e_0=\phi_0-s\phi_1+\phi_2,\qquad e_1=\phi_1-s\phi_2+\phi_3,\qquad
u_\pm=\frac{e_0\pm e_1}{\|e_0\pm e_1\|_2}.
\]
They span \(E=\operatorname{span}_{\mathbb C}\{u_+,u_-\}\) with \(q(u_\pm,u_\pm)=\mu_\pm>7\) (T1C), and \(q(u_+,u_-)=0\). Pairing \(q(f,g)=W(f*\widetilde g)\) is linear in the first argument and uses the **original complete signed Weil form**. The physical \(L^2(I)\) inner product is the denominator throughout.

T1D-J defined the reflected, shifted three-cell unit vectors
\[
y_{2,R}(x)=\frac{e_0(x-\log2)}{\sqrt{JN_y}},\qquad
y_{2,L}(x)=y_{2,R}(-x),\qquad N_y=2+s^2=22/3.
\]
Add, **without changing the bump shape**, a second translated reflection pair:
\[
\boxed{y_{4,R}(x)=\frac{e_0(x-\log4)}{\sqrt{JN_y}},\qquad y_{4,L}(x)=y_{4,R}(-x).}\tag{K.1}
\]
The furthest cell lies at \(\tfrac12\log3+\log4<2\), rigorously from \(4\sqrt3<e^2\), and the same holds after reflection. All 16 individual source cells (four original, three each from the four translated probes) are **pairwise disjoint**: their minimum multiplicative center ratio is \(9/8\), and \(\ell<\log(9/8)\). Every source is smooth, compactly supported in \(I\), and satisfies **both** Mellin pole moment constraints exactly, because translation multiplies each moment by a nonzero factor while reflection exchanges the two moments. Hence the six vectors are admissible, and all pole terms in their pairwise Weil form vanish by exact annihilation.

Define \(v_{h,\pm}=(y_{h,R}\pm y_{h,L})/\sqrt2\) for \(h=2,4\). The six vectors \(u_+,u_-,v_{2,+},v_{2,-},v_{4,+},v_{4,-}\) form a **physical orthonormal** family. Let
\[
\boxed{V_{2,4}=\operatorname{span}_{\mathbb C}\{v_{2,+},v_{2,-},v_{4,+},v_{4,-}\}\subset F,
\quad\mathscr E_6=E\oplus V_{2,4},\ \dim\mathscr E_6=6,}\tag{K.2}
\]
where \(F=C_c^\infty(I)\cap(E\oplus P)^\perp\), \(P=\operatorname{span}\{e^{x/2}\mathbf1_I,e^{-x/2}\mathbf1_I\}\). The smooth sources lie in \(F\); the nonsmooth windowed pole functions are **not** used as Weil tests.

## 2. Exhaustive original prime-power arithmetic

For the fixed support radius \(a=2\), the exact cutoff is \(n<e^4\), i.e. the 24 genuine prime powers \(2\le n\le54\). Represent each source center as \(x=(a/2)\log3+b\log2\) with \((a,b)\in\mathbb Z^2\); all relevant center ratios \(3^{(a-a')/2}2^{b-b'}\) are rational because every difference \(a-a'\) is even. Checking the strict finite-width windows using \(e^\ell\le(1-\ell)^{-1}\) gives exactly:

| Pair of source groups | All cross-group prime-power coincidences |
|---|---|
| Original base cells and shifted \(2\)-cells | \(2\) (three equal-index pairs) |
| Original base cells and shifted \(4\)-cells | \(4\) (three equal-index pairs) |
| Right \(2\)-cells and reflected left \(2\)-cells | \(4\) (two pairs) |
| Right \(4\)-cells and reflected left \(4\)-cells | \(16=2^4\) (two pairs) |
| Right \(2\)-cells and right \(4\)-cells | \(2\) (three pairs) |
| Right \(2\)-cells and reflected left \(4\)-cells | \(8=2^3\) (two pairs) |
| Reflected left \(2\)-cells and right \(4\)-cells | \(8=2^3\) (two pairs) |
| Reflected left \(2\)-cells and reflected left \(4\)-cells | \(2\) (three pairs) |

Within each translated three-cell block, the only off-diagonal prime powers remain \(3\) and \(9\). **No composite with \(\Lambda(n)=0\) has been counted as a prime.** All source prime-power collisions are enumerated directly by `src/audits/canonical_return_t1d_k_extended_probes.py`, which also constructs the entire finite six-mode prime Gram algebraically from these 16 cells.

Write \(\kappa_n=\Lambda(n)/\sqrt n\), \(D_+=2-2s+s^2\), \(D_-=2+2s+s^2\), and \(N_\pm=(44\mp16\sqrt3)/3\). Exact parity reduces the complete six-dimensional Weil Gram to an orthogonal sum of **two \(3\times3\) real Hermitian blocks**, in basis \((u_\sigma,v_{2,\sigma},v_{4,\sigma})\), \(\sigma=\pm\). Their genuine **prime-only** off-diagonal entries are
\[
\boxed{\begin{aligned}
q_{\rm prime}(u_\sigma,v_{2,\sigma})&=-\frac{\sqrt2\,\kappa_2D_\sigma}{\sqrt{N_\sigma N_y}},\\
q_{\rm prime}(u_\sigma,v_{4,\sigma})&=-\frac{\sqrt2\,\kappa_4D_\sigma}{\sqrt{N_\sigma N_y}},\\
q_{\rm prime}(v_{2,\sigma},v_{4,\sigma})&=-\kappa_2+\epsilon_\sigma\frac{2s}{N_y}\kappa_8,
\quad \epsilon_+=1,\ \epsilon_-=-1.
\end{aligned}}\tag{K.3}
\]
The translated parity-probe diagonals each have the original three-cell prime contribution \(7\log3/11\) and a further reflected interaction \(\epsilon_\sigma(2s/N_y)\kappa_4\) for \(h=2\), or \(\epsilon_\sigma(2s/N_y)\kappa_{16}\) for \(h=4\). The **Gamma terms are retained separately**; (K.3) is never substituted for the full Gram.

## 3. Complete Gamma control and six-dimensional positivity

For disjoint sources, Suzuki's archimedean cross-kernel is \(-w(|x-y|)\), \(w(t)=e^{-t/2}/(1-e^{-2t})\). All cross-group cell separations exceed \(\log(9/8)-\ell\), so
\[
w(t)\le\frac1{1-(8/9)^2/(1-2\ell)}<5.
\]
The base-to-2 translations have stronger minimum center ratio \(3/2\), giving \(w<2\). Because \((\int|\phi|)^2/J\le\ell\), the Gamma mixed pairing between two coefficient combinations is bounded by \(5\ell\) times the product of their coefficient \(\ell^1\) masses divided by their physical normalization factors. For the two translated three-cell parity diagonals, J's earlier estimate still applies: \(q(v_{h,\sigma},v_{h,\sigma})>36/5\) for **both** \(h=2\) and \(h=4\), because the \(h=4\) reflected prime coupling is smaller (\(\kappa_{16}<\kappa_4\)) and its Gamma cross-kernel is no larger than the safe \(3\ell(2+s)^2/N_y<1/1000\) used for \(h=2\). The original T1C direct Gamma diagonal certificate \(g_0/J>679/100\) is imported as a source-level premise, not silently re-proved here.

Exact-rational outward enclosures for the **prime magnitudes** and **Gamma absolute errors** in the four selected/probe scalar couplings are:

| Source pair | Exact prime magnitude interval | Gamma absolute upper |
|---|---:|---:|
| \(u_+,y_{2,R}\) | \([0.21085901,0.21085902]\) | \(<0.00063091\) |
| \(u_+,y_{4,R}\) | \([0.14909983,0.14909984]\) | \(<0.00157726\) |
| \(u_-,y_{2,R}\) | \([0.44245336,0.44245337]\) | \(<0.00056106\) |
| \(u_-,y_{4,R}\) | \([0.31286177,0.31286178]\) | \(<0.00140264\) |

Thus the **complete** (prime + Gamma) mixed terms satisfy
\[
\boxed{|q(u_+,v_{2,+})|<\frac{\sqrt2}{4},\quad
|q(u_+,v_{4,+})|<\frac{\sqrt2}{5},\quad
|q(u_-,v_{2,-})|<\frac{\sqrt2}{2},\quad
|q(u_-,v_{4,-})|<\frac{2\sqrt2}{5}.}\tag{K.4}
\]
The \(v_{2,\sigma}\)-to-\(v_{4,\sigma}\) prime expression (K.3), together with its two signed Gamma cross-correlations (each \(\le5\ell(2+s)^2/N_y\)), gives
\[
\boxed{|q(v_{2,\sigma},v_{4,\sigma})|<\frac23.}\tag{K.5}
\]
By strict Gershgorin/diagonal dominance in each of the two physical-orthonormal \(3\times3\) blocks, with diagonals \(>7,>36/5,>36/5\), and \(\sqrt2<3/2\), equations (K.4)–(K.5) imply
\[
\boxed{\forall\,0\ne f\in\mathscr E_6:\qquad q(f,f)>\frac{11}{2}\|f\|_2^2.}\tag{K.6}
\]
This is a genuinely **six-dimensional** compact smooth *two-Mellin-pole-null* positivity result for the complete original signed Weil form, not a full-form positivity theorem. It retains \(\Lambda(2),\Lambda(4),\Lambda(8),\Lambda(16)\) with the correct signs and weights, in addition to \(3,9,27\) and the complete Gamma distribution. These source tests are not known to be Connes's conditioned Sonin/prolate vectors.

## 4. The enlarged four-dimensional complement probe passes its Schur inequality

Write \(\mathcal F=(P\oplus E)^\perp\subset L^2(I)\), \(F=\mathcal D^0\cap\mathcal F\), and \(r_\sigma=\Pi_\mathcal F\mathcal K_2u_\sigma\). The exact positive rank-two penalty is \(\mathcal R=\sum_\sigma|r_\sigma\rangle\langle r_\sigma|/\mu_\sigma\), with \(\mu_\sigma>7\). Because \(V_{2,4}\subset F\) is physically orthonormal and parity-diagonal,
\[
\boxed{
M_\sigma:=\|\Pi_{V_{2,4}}r_\sigma\|_2^2
=2\bigl(|q(u_\sigma,y_{2,R})|^2+|q(u_\sigma,y_{4,R})|^2\bigr).
}\tag{K.7}
\]
The full Weil mixed Gamma error bounds above yield **exact-fraction, outward-propagated** captured-residual intervals
\[
\boxed{\begin{aligned}
0.13191753&\le M_+\le0.13486315,\\
0.58455122&\le M_-\le0.59004781.
\end{aligned}}\tag{K.8}
\]
The convenient rigorous rational inequalities are
\[
\boxed{\frac{13}{100}<M_+<\frac7{50},\qquad
\frac{29}{50}<M_-<\frac35.}\tag{K.9}
\]
Consequently
\[
\boxed{\|\Pi_{V_{2,4}}\mathcal R\Pi_{V_{2,4}}\|
=\max\{M_+/\mu_+,M_-/\mu_-\}<\frac3{35}.}\tag{K.10}
\]
Also, both \(2\times2\) parity blocks of \(q|_{V_{2,4}}\) have diagonal \(>36/5\) and off-diagonal \(<2/3\), hence \(q|_{V_{2,4}}>13/2\) in the physical metric. Therefore the **true source Schur remainder is positive on all four complementary probes**:
\[
\boxed{
q(y,y)-\langle\mathcal R y,y\rangle>\frac{32}{5}\|y\|_2^2
\quad\forall\,0\ne y\in V_{2,4}.}\tag{K.11}
\]
Again, this is a test on a **four-dimensional subspace of** \(F\), not on all of \(F\).

## 5. An explicit lower bound for how much complement coupling remains unmeasured

Let \(r_\sigma^{\perp}= (I-\Pi_{V_{2,4}})r_\sigma\). Because \(V_{2,4}\subset\mathcal F\), the source **identity**
\[
\boxed{\|r_\sigma^\perp\|_2^2=\beta_\sigma^2-M_\sigma}\tag{K.12}
\]
is exact and uses only physical orthogonality. To turn it into a numerical lower bound, **we separately import** the analytic two-sided residual estimates from T1D-I, \(\beta_+^2>16/5\), \(\beta_-^2>12/5\). Together with (K.9) they yield
\[
\boxed{\|r_+^\perp\|_2^2>\frac{153}{50}=3.06,\qquad
\|r_-^\perp\|_2^2>\frac95=1.8.}\tag{K.13}
\]
The unmeasured mixed norm remains substantial even after including both prime-2 and prime-4 translated reflection pairs. This conclusion is **conditional on the correctness of T1D-I's infinite-dimensional analytic estimates**, not independently certified by the new finite symbolic script. Even without importing T1D-I, (K.7)–(K.11) remain valid finite-source results.

The same separation has a convenient operator-norm formulation. T1D-I reported \(\|\mathcal R\|>320/1393\); (K.10) gives
\[
\boxed{
\|\mathcal R-\Pi_{V_{2,4}}\mathcal R\Pi_{V_{2,4}}\|
\ge \|\mathcal R\|-\|\Pi_{V_{2,4}}\mathcal R\Pi_{V_{2,4}}\|
>\frac{1003}{6965}\approx0.1440.}\tag{K.14}
\]
This does **not** say the complement form has a negative direction. It says that the strong source-defined mixed interactions cannot be fully captured by these particular natural arithmetic probes.

## 6. Transfer boundary, falsifiers and next research gate

**PASS on the smooth original Weil source (with separately cited analytic inputs):** exact pole-null/six-direction physical orthonormality, complete prime-power census and symbolic 6×6 prime Gram, rigorous Gamma cross estimates, complete signed Weil positivity \(q|_{\mathscr E_6}>11/2\), and source Schur domination on the four-dimensional authentic complement probe \(V_{2,4}\). The independent checker evaluates exact rational inequalities and full finite prime matrix (not every smooth integral). The quantitative *outside* residual bound (K.13) additionally invokes T1D-I, as clearly labeled.

**Falsify/review at:** actual support margins for the reflected \(4\)-translated cells, the \(8\)- and \(16\)-prime interaction masks, negative von Mangoldt signs, the fixed \(\mu_\pm>7\) input, the global Gamma cross-weight majorant, and T1D-I's non-compact logarithmic Fourier-tail estimate. No original Weil-Lab/RH repository is altered.

**OPEN:** sign or coercivity of \(q\) on *all* \(F\), sharp two-sided \(\beta_\pm\) quadrature, basis-complete probe estimates, Connes's H2′/R1/R2 conditioned Sonin space and operator-domain correspondence, two source pole directions in the full Weil test criterion, global Weil positivity, and RH. No spectral-level inference from a finite positive compression.

**Next CR-INT-T1D-L:** construct and certify a **frequency-conditioned infinite-dimensional complement** rather than adding a third isolated prime translate. Determine whether a projected Sonin/prolate condition defines a comparable source domain under both metrics, or document a strict incompatibility. If numerical tests are used, distinguish finite-matrix spectral diagnostics from actual complement coercivity.
