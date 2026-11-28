# CR-INT-T1D-J — Genuine two-dimensional complement probe, prime-4 carrier, and four-dimensional Weil positivity

**Canonical Return · independent signed-Weil smooth-core research · 2026-11-23 (America/Los_Angeles).**  
**Status:** **DERIVED / ANALYTIC**, with direct prime-power collision enumeration, exact SymPy compressed-matrix assertions, and outward-rational conservative inequalities **AUDIT-PASS**. This does **not** establish coercivity of the *entire* infinite-dimensional complement, a spectral theorem, a completed-domain Suzuki–Connes equivalence, a Lean certificate, or RH. It follows [T1C](../archive/space/10-t1c.md), [T1D-D](CR-INT-T1D-D.md), [T1D-E](CR-INT-T1D-E.md), [T1D-H](CR-INT-T1D-H.md), and [T1D-I](CR-INT-T1D-I.md).

## 1. Genuine source and exact geometry

Use the complete **signed original Weil form** from [Suzuki 2026, arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4), with the two Mellin-pole conditions. Fix the physical interval \(I=(-2,2)\), the exact **even, nonnegative, smooth** triangle–mollifier \(\phi\) from T1C (width \(\ell=10^{-4}\)), its mass \(J=\|\phi\|_2^2>0\), \(d=\log3\), and \(s=4/\sqrt3\). Let \(\phi_j\) be its translations to centers \(x_j=(j-3/2)d\), \(0\le j\le3\), and put

\[
 e_0=\phi_0-s\phi_1+\phi_2,\qquad
 e_1=\phi_1-s\phi_2+\phi_3,\qquad
 u_\pm=\frac{e_0\pm e_1}{\|e_0\pm e_1\|_2}.
\]

As proved and audited in T1B/C, \(E=\mathrm{span}_{\mathbb C}\{u_+,u_-\}\) is pole-null and physically orthonormal, with complete Weil values \(q(u_\pm,u_\pm)=\mu_\pm>7\) (T1C's stronger numerical intervals are not required below). With the pairing linear in its first argument, write \(q(f,g)=W(f*\widetilde g)\).

The key change is that we now test **two independent genuine smooth complementary functions**, rather than merely naming a general \(y\in F\). Define

\[
 \boxed{y_R(x)=\frac{e_0(x-\log2)}{\sqrt{JN_y}},\qquad
 y_L(x)=y_R(-x)=\frac{e_1(x+\log2)}{\sqrt{JN_y}},\qquad
 N_y=2+s^2=\frac{22}{3}.}\tag{J.1}
\]

The three centers of \(y_R\) are \(x_0+\log2,x_1+\log2,x_2+\log2\), and those of \(y_L\) are their reflections. All six lie strictly inside \(I\). Their supports are pairwise disjoint (and disjoint from all four original cells) at width \(10^{-4}\), so \(\|y_R\|_2=\|y_L\|_2=1\) and \(\langle y_R,y_L\rangle=0\). Both have **both Mellin pole moments exactly zero**, since the three-cell annihilator \(q(z)=1-sz+z^2\) is preserved by translation and reflection. Therefore

\[
 \boxed{V=\mathrm{span}_{\mathbb C}\{y_R,y_L\}\subset F:=\mathcal D_2^0\cap E^\perp,\qquad V\perp E.}\tag{J.2}
\]

This source identity is exact; the pole functions \(e^{\pm x/2}\mathbf1_I\) are used for the *physical* complement only and are not substituted as smooth Weil tests.

## 2. Exhaustive prime-power census: new prime 4, not an assumed approximation

For cells with center separation \(\log r\), the prime \(n\) term can occur only when \(|\log(r/n)|<\ell\) (with \(r\) replaced by \(r^{-1}\) if necessary). The exact rational bound \(e^\ell\le(1-\ell)^{-1}\), combined with enumeration of **all 24 genuine prime powers \(n\le54\)**, establishes these *complete* cross-correlation masks:

| Pair | Only activated prime-power coincidences | How many |
|---|---|---:|
| Original 4-cell source against \(y_R\) | \(2\), matching the first three original cells | 3 |
| Original 4-cell source against \(y_L\) | \(2\), matching the last three original cells | 3 |
| \(y_R\) against itself (off-diagonal) | \(3\) between neighbors, \(9\) between extremes | 4 and 2 oriented entries |
| \(y_R\) against \(y_L\) | **\(4=2^2\)**, in precisely two ordered cell pairs | 2 |

For the last row, the ratio of a right-cell center to a reflected left-cell center is \(4\cdot3^{j+k-3}\), \(j,k=0,1,2\); the five possible ratios are \(4/27,4/9,4/3,4,12\). Only \(4\) is a genuine prime power in the corresponding narrow window. **In particular, the other potential rational ratios do not secretly introduce composite von Mangoldt terms.** The audit constructs the prime-only \(4\times4\) Gram matrix directly from these collision masks and checks its exact symbolic entries.

Set \(\kappa_n=\Lambda(n)/\sqrt n\), so \(\kappa_2=\log2/\sqrt2\), \(\kappa_4=\log2/2\), \(\kappa_3=\log3/\sqrt3\), and \(\kappa_9=\log3/3\). Write \(D_+=2-2s+s^2\), \(D_-=2+2s+s^2\), and \(N_\pm=(44\mp16\sqrt3)/3\). The **exact genuine-prime pieces** of the normalized Gram are

\[
\boxed{\begin{aligned}
 q_{\rm prime}(u_\pm,y_R)&=-\frac{\kappa_2D_\pm}{\sqrt{N_\pm N_y}},\\
 q_{\rm prime}(u_+,y_L)&=q_{\rm prime}(u_+,y_R),\\
 q_{\rm prime}(u_-,y_L)&=-q_{\rm prime}(u_-,y_R),\\
 q_{\rm prime}(y_R,y_R)&=\frac{4s\kappa_3-2\kappa_9}{N_y}=\frac7{11}\log3,\\
 q_{\rm prime}(y_R,y_L)&=\frac{2s\kappa_4}{N_y}>0.
\end{aligned}}\tag{J.3}
\]

The complete Weil form also includes **the full Gamma diagonal and every Gamma cross term**. The pole terms vanish by exact source nullity, not by dropping them. Equation (J.3) alone is *not* called the full Weil Gram.

## 3. Complete-form, non-asymptotic lower bounds for V

Let \(g_0\) be T1C's full single-bump Gamma diagonal and \(A_k\) its genuine cross-Gamma correlations. The exact diagonal Weil value on either complementary vector is

\[
\boxed{q(y_R,y_R)=q(y_L,y_L)=\frac{g_0}{J}+\frac7{11}\log3+\frac{4s(A_1/J)-2(A_2/J)}{N_y}.}\tag{J.4}
\]

Here \(A_k\) depends on the frozen smooth \(\phi\). Even without using nonnegativity, \(|A_k|/J<\ell\): Cauchy–Schwarz gives \((\int|\phi|)^2/J\le\ell\), and the Gamma weight \(w(t)=e^{-t/2}/(1-e^{-2t})<1\) for all cross-support separations \(t\ge\log3-\ell>\log2\). The independently executed exact-fraction T1C audit gives \(g_0/J>679/100\); an atanh-series check gives \(\log3>109/100\). Since \(2<s<7/3\), \(N_y>7\), and \(\ell=10^{-4}\), this yields

\[
\boxed{q(y_R,y_R)>\frac{679}{100}+\frac7{11}\frac{109}{100}-\frac1{5000}>\frac{149}{20}=7.45.}\tag{J.5}
\]

For the mixed complement pair, every center separation is at least \(\log(4/3)\); after accounting for the bump width it exceeds \(\log(5/4)\), where the monotone Gamma kernel satisfies \(w(t)<3\). Normalized Cauchy–Schwarz therefore gives

\[
 \bigl|q_\Gamma(y_R,y_L)\bigr|
 <\frac{3\ell(2+s)^2}{N_y}<\frac1{1000}.
\]

The prime-4 part obeys \((2s/N_y)\kappa_4< (2/3)(7/20)=7/30\). Hence

\[
\boxed{|q(y_R,y_L)|<\frac7{30}+\frac1{1000}<\frac14.}\tag{J.6}
\]

The reflected even and odd normalized probes

\[
\boxed{v_+=\frac{y_R+y_L}{\sqrt2},\qquad v_-=\frac{y_R-y_L}{\sqrt2}}\tag{J.7}
\]

form a *physical orthonormal basis of V*. The source is reflection invariant, so \(q(v_+,v_-)=0\), and (J.5)–(J.6) prove the **genuine complementary compression margin**

\[
\boxed{q(v_\pm,v_\pm)>7.45-0.25=\frac{36}{5}=7.2,\qquad
q(y,y)>\frac{36}{5}\|y\|_2^2\quad(0\ne y\in V).}\tag{J.8}
\]

This is a two-dimensional **subspace of the authentic infinite-dimensional complement**, not a positivity theorem for that whole complement.

## 4. Exact cross-block test and a new 4D positive compression

All selected–complement pairings use the **original full signed Weil form**. The closest separated cell centers are \(\log(3/2)\) apart, so their Gamma cross weight is \(<2\). With \(\int|\phi|\le\sqrt{\ell J}\), the complete cross-Gamma magnitude is bounded by

\[
 |q_\Gamma(u_\sigma,y_R)|
 <\frac{2\ell\,\|c_\sigma\|_1(2+s)}{\sqrt{N_\sigma N_y}}<\frac1{1000}.
\]

For the exact prime-2 term in (J.3), use \(\log2/\sqrt2<1/2\), \(s<7/3\), \(\sqrt{N_+N_y}>6\), and \(\sqrt{N_-N_y}>13\). The independently checked inequalities are

\[
\boxed{|q(u_+,y_R)|<\frac14,\qquad |q(u_-,y_R)|<\frac12.}\tag{J.9}
\]

By parity, \(q(u_+,v_-)=q(u_-,v_+)=0\) **exactly**, while

\[
 q(u_+,v_+)=\sqrt2\,q(u_+,y_R),\qquad
 q(u_-,v_-)=\sqrt2\,q(u_-,y_R).
\]

Thus the **complete signed Weil form on the four-dimensional pole-null source \(E\oplus V\)** is the orthogonal sum of two real Hermitian \(2\times2\) blocks, each with selected diagonal \(\mu_\sigma>7\), complementary diagonal \(>36/5\), and mixed entry of magnitude \(<\sqrt2/2<3/4\). The elementary \(2\times2\) diagonal-dominance inequality gives the *new* strictly positive four-dimensional certificate

\[
\boxed{\forall f\in E\oplus V\setminus\{0\},\qquad q(f,f)>6\|f\|_2^2.}\tag{J.10}
\]

This is **not** a statement about the other infinitely many source directions; it adds exactly two legitimate physical-complement probes while retaining all poles, genuine prime powers, and Gamma contributions.

## 5. Schur diagnostic: positivity on V does not close the complement

Recall T1D-E's positive rank-two Schur penalty \(\mathcal R=\sum_\sigma |r_\sigma\rangle\langle r_\sigma|/\mu_\sigma\) on \(\mathcal F=(E\oplus P)^\perp\). Its compression onto V is diagonal in the parity basis, and (J.9) with \(\mu_\pm>7\) yields

\[
\boxed{0\le\|\Pi_V\mathcal R\Pi_V\|<\frac{2(1/2)^2}{7}=\frac1{14}.}\tag{J.11}
\]

Combining this with (J.8) gives a **strictly positive Schur remainder on V**:

\[
\boxed{q(y,y)-\langle\mathcal R y,y\rangle>7\|y\|_2^2\qquad(0\ne y\in V).}\tag{J.12}
\]

This is a real, finite-dimensional *test* of the exact infinite-complement inequality, not a proof of that inequality for every \(y\in F\).

There is an additional diagnostic **conditional on the independent T1D-I two-sided residual analysis**: its reported full-complement lower bound \(\|\mathcal R\|>320/1393\) is strictly larger than the present \(\|\Pi_V\mathcal R\Pi_V\|<1/14\). If the T1D-I interval is upheld under external review, the explicit prime-2 translated probe space is **not** a norm-maximizing subspace of the true residual penalty, with the quantitative norm separation

\[
\|\mathcal R-\Pi_V\mathcal R\Pi_V\|>\frac{320}{1393}-\frac1{14}=\frac{63}{398}.
\]

This last comparison *depends* on T1D-I and is not independently re-certified here; (J.10)–(J.12) do not depend on it.

## 6. Transfer / stopping conditions / verification

**New analytic source result:** two explicit, linearly independent smooth \(y_R,y_L\in F\) carry exactly prime-2 selected cross-couplings and a genuine prime-4 internal interaction; their complete Weil compression obeys (J.8), their source-defined Schur remainder obeys (J.12), and the extended **four-dimensional** pole-null span satisfies (J.10). All statements are at the original **fixed T1C width \(10^{-4}\)** and physical interval \((-2,2)\).

**Machine verification:** `src/audits/canonical_return_t1d_j_complement_probe.py` re-executes the original T1C exact rational enclosure, independently enumerates all 24 prime powers and the cross-support masks, constructs all prime Gram entries by direct collision enumeration, asserts the exact Mellin null matrix and native norm factors, and proves the conservative rational margin chain. Gamma integral representation, support-source passage, and parity are written analytical arguments, **not Lean compilation**. These statuses must not be conflated.

**Unresolved:** a positive bound on \(q|_F\) for **all** infinitely many directions, a bound on the operator Schur penalty relative to the actual form (rather than on V), Connes's conditioned Sonin/prolate subspace and H2′/R1/R2 hypotheses, completed-domain transport to Suzuki/Connes, the two remaining pole directions of the full Weil criterion, and RH. The newly verified four-dimensional positive compression **cannot be iterated to prove global positivity without a uniform and complete basis/limit theorem**.

**Next — CR-INT-T1D-K:** compute the physical residual *outside* the current complement probes, e.g. the orthogonal residual \((I-\Pi_V)r_\sigma\), and test a rigorously specified basis of \(\mathcal F\) or a genuine frequency-conditioned complement. A successful finite probe does not license extrapolation to full \(F\); investigate explicit spectral/Schur estimates on the untested orthogonal directions before any Suzuki–Connes equivalence conjecture is promoted.