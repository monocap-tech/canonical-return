# CR-ALP-019 — Four-prime smooth pole-null family: complete original Weil 4×4 Schur certificate

**Canonical-Return | 2026-11-23 (America/Los_Angeles) | research-only / draft.**

**Standing:** exact source-contract comparison; analytical finite-subspace theorem with independently executed **pure-Fraction rational audit PASS**. **NOT** new Lean compilation, F41 numerical H/V entry certification, positivity on the entire old support core, arbitrary-support extension, Connes/Suzuki completed-domain admission, or RH.

## 0. Independently reopened source custody and a necessary type distinction

The receiving **actual original** IV010 and IV011 definitions were reopened from `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`:

- `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportExtension.lean`, Git blob `5ef51fd37337c466620eafccee2a27b77955746e`: `neutralIV010Mixed` is the literal `Zeta23.EF.literatureRHS` of the physical first-slot-conjugated correlation. The exact complex Hermitian old/new block and its **conditional** Schur theorem appear in `neutralIV010OldPlusNew_realBlock` and `neutralIV010OriginalOneNewSource_nonnegative_of_crossBound`.
- `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`: `neutralIV011MixedCorrelation_eq_schwartzIntegral`, `neutralIV011OriginalMixedPrime_eq_finite`, `neutralIV011OriginalMixed_arithmeticSplit` and `neutralIV011OriginalOneNew_nonnegative_iff_crossBound` preserve exactly *pole − both original prime shifts + complete archimedean Gamma*, with the **mixed** natural cutoff `log n > a_old + a_new`.
- **Separately**, the historical original physical F41 source was fetched by its actual immutable commit `monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5`. File `WeilDefect/Arithmetic/ActualZetaOriginalPhysicalSolveResidual.lean`, blob `2cf06688dc7e63bbb0c797102260559884325855`, explicitly declares
  `neutralNF606PhysicalH(i,j) = ∫_{-11/10}^{11/10} nf545OriginalSelectedRealGenerator(i,x) * nf545OriginalSelectedRealGenerator(j,x) dx`,
  and `neutralNF606PhysicalV(i,j) = ∫ nf545OriginalSelectedRealGenerator(i,x) * neutralNF587CriticalGenerator(j,x) dx`.

**Type warning:** the F41 historical **H is a native physical L² Gram**, with 242 old generators, and its V is 242×4 native physical mixed data; it is **not** a signed Weil 242×242 matrix. In this note **M** denotes the native physical Gram and **H_W** denotes the *signed Weil quadratic form* on **newly constructed packet directions**. Neither matches F41's frozen generator inventory without an actual source-defined embedding and complete mixed-moment computation. Prior CR notes should not be read as identifying those matrices.

The original Weil Gamma/contact dictionary and the physical source convention were reopened in `IV_CR_TRANSFER_II_ORIGINAL_RADIUS_OPERATOR.md` (Phase IV analysis, 23 November 2026); the complete signed distribution is independently retained here, with no substitute positive form.

## 1. Four explicitly centered, smooth, two-pole-null sources

Fix

\[
\boxed{b=10^{-8},\quad \eta=10^{-14},\quad w=2(b+\eta).} \tag{19.1}
\]

Take **any** real even nonnegative mass-one smooth mollifier `ρ_η∈C_c^∞(-η,η)` and set

\[
\begin{aligned}
 g_b(u)&=e^{u/2}{\bf1}_{[0,b]}(u)-e^{u/2-b}{\bf1}_{[b,2b]}(u),\\
 g&=\rho_\eta*g_b,\qquad G=\|g\|_2^2.
\end{aligned}
\]

For each distinct prime `p∈P={2,3,5,7}`, put

\[
\boxed{h_p(u)=g(u+\log p+b)-p\,g(u-\log p+b).} \tag{19.2}
\]

The two packets lie close to `−log p` and `+log p`, each with width at most `w`. The original two Mellin moments obey `F_+(g_b)=0`, and for the uncentered pair `g_b−p T_{2log p}g_b`, `F_−` acquires the factor `1−p e^{−log p}=0`. Convolution and common translation preserve exact zero moments. Therefore

\[
\boxed{F_+(h_p)=F_-(h_p)=0,\qquad h_p\in C_c^\infty(\mathbb R),\quad h_p\ne0.} \tag{19.3}
\]

Each is a Schwartz source. Its centered support radius is `a_p=log p+b+η`. Exact rational log enclosures certify the strict inclusion chain:

| Prime | Certified physical radius containing the source |
|---|---|
| `p=2` | `a_2<4/5` |
| `p=3` | `a_3<11/10` — the reviewed original Phase IV radius |
| `p=5` | `a_5<17/10` |
| `p=7` | `a_7<2` |

Thus `h_3`, `h_5`, `h_7` are genuinely new with respect to the preceding radius, rather than merely rescalings of an old function. The eight individual packets are mutually disjoint, so

\[
\boxed{M_{pq}=\langle h_p,h_q\rangle_{L^2}
=\delta_{pq}(1+p^2)G.}\tag{19.4}
\]

On the positive interior `[η,b−η]`, `g(u)=m_ηe^{u/2}` with `m_η=∫ρ_η(v)e^{−v/2}dv≥1` (evenness and Jensen/cosh). Therefore

\[
\boxed{G\ge b-2\eta>0.}\tag{19.5}
\]

For the chosen prime order the diagonal physical masses are `(5G,10G,26G,50G)`; **this is the native source metric** and not a frame renormalization.

## 2. Entire original pole and prime pieces, including mixed terms

Every crossed two-pole form vanishes since each source kills both Mellin moments. For a single `p`, the autocorrelation `k_p=h_p*\widetilde h_p` is supported in the three `w`-neighborhoods of `0, ±2log p`, and

\[
\boxed{k_p(2\log p)=-pG.}\tag{19.6}
\]

The exact rational support audit shows only the **prime power `n=p²`** survives the original positive-frequency sample. Hence the full original signed prime term on that diagonal is

\[
\boxed{\mathcal W_{\rm prime}(h_p,h_p)
=-2\frac{\log p}{p}\,(-pG)=2(\log p)G>0.}\tag{19.7}
\]

For two **distinct primes** `p<q`, the authentic first-slot-conjugated cross correlation is supported only near

\[
\boxed{\pm\log(q/p),\qquad \pm\log(pq).}\tag{19.8}
\]

The ratio `q/p` is never an integer, and `pq` is **not a prime power** (its von Mangoldt value is zero). The width `w` is smaller than every distance from these peaks to a genuine prime-power log. The auditor verifies **every** prime power `n≤50` using the exact inequality

\[
|\log r-\log n|\ge\frac{|r-n|}{\max(r,n)}
\]

and independently excludes `n≥51`, as all peaks lie at most at `log 49`. Thus **both signed frequency contributions are accounted for**:

\[
\boxed{\mathcal W_{\rm prime}(h_p,h_q)=0\quad(p\ne q).}\tag{19.9}
\]

In particular, this is not a claim that the prime side is positive **termwise** for an arbitrary compact-supported source; it is an exact geometric cancellation for these particular packets.

## 3. Complete original Gamma/contact term gives a four-dimensional quantitative margin

The full Gamma distribution on an even real autocorrelation is, with the literal original convention,

\[
\mathcal W_\Gamma(k)
=m_\Gamma(0)k(0)+2\int_0^\infty[k(0)-k(t)]w_\Gamma(t)\,dt,\quad
w_\Gamma(t)=\frac{e^{-t/2}}{1-e^{-2t}},
\]

and `m_Γ(0)=−γ−3log2−π/2−logπ>−6`. The exact rational auditor establishes the elementary numerical brackets used below, including a Machin rational interval for `π`.

**Diagonal.** By Cauchy–Schwarz, `k_p(t)≤k_p(0)=(1+p²)G`; the difference is everywhere nonnegative. The packet gap implies `k_p(t)=0` for `t∈[w,1]`, since `2log2=log4>1+w`. On that interval `w_Γ(t)≥e^{−1/2}/(2t)`. The exact rational bounds

\[
e^{-1/2}>3/5,\qquad\log(1/w)>16,\qquad m_\Gamma(0)>-6
\]

give

\[
\boxed{\mathcal W_\Gamma(h_p,h_p)>
\frac{18}{5}(1+p^2)G.}\tag{19.10}
\]

The original prime contribution (19.7) is additionally positive, so every **complete signed** diagonal is strictly above this same threshold.

**Mixed.** Since `h_p,h_q` are disjoint, their mixed correlation at zero vanishes and the original Gamma contact term drops out. The positive mixed-lag peaks satisfy

\[
 t\ge\log(7/5)-w>1/4.
\]

For `t>1/4`, `w_Γ(t)<3` because `e^{1/2}>3/2`. The whole crossed original Gamma integral is therefore bounded via convolution Young's inequality:

\[
\begin{aligned}
|\mathcal W_\Gamma(h_p,h_q)|
&<3\|h_p\|_1\|h_q\|_1\\
&<27(1+p)(1+q)b^2\\
&\le1296b^2.
\end{aligned}\tag{19.11}
\]

Here `||g||₁≤||g_b||₁=4sinh(b/2)<3b`; every pole and mixed prime term is *exactly zero*, so (19.11) is a bound on the **entire original signed mixed Weil entry** `B_pq`.

Let `H_W=(W(h_p,h_q))_{p,q∈P}`, in the correct first-slot-antilinear Hermitian convention. Put `\widetilde H=M^{-1/2}H_WM^{-1/2}`. Since the two smallest distinct masses are `5G` and `10G`,

\[
\frac{|B_{pq}|}{\sqrt{M_{pp}M_{qq}}}
<\frac{1296b^2}{7G}
\le\frac{1296b^2}{7(b-2\eta)}.
\]

Every row has at most three off-diagonal entries. The pure-Fraction audit proves

\[
\boxed{\sum_{q\ne p}|\widetilde H_{pq}|
<\frac{3888b^2}{7(b-2\eta)}<10^{-4}.}\tag{19.12}
\]

The diagonal entries exceed `18/5`. Hermitian Gershgorin—or direct expansion with `2|z_i z_j|≤|z_i|²+|z_j|²`—therefore gives

\[
\boxed{
H_W>\left(\frac{18}{5}-\frac1{10000}\right)M
=\frac{35999}{10000}M
\quad\text{on }\mathbb C^4.}
\tag{19.13}
\]

Equivalently, for every nonzero complex vector `z=(z_2,z_3,z_5,z_7)`,

\[
\boxed{
\mathcal W\!\left((\sum_p z_ph_p)*\widetilde{(\sum_p z_ph_p)}\right)
>\frac{35999}{10000}\Big\|\sum_p z_ph_p\Big\|_2^2>0.}
\tag{19.14}
\]

This is **genuine, complete, original signed Weil** positivity on a **selected four-dimensional smooth pole-null physical subspace**. It does **not** prove positivity on the containing infinite-dimensional window.

### 3a. Three actual successive finite Schur extension certificates

Let `E_j=span{h_2,…,h_{p_j}}` in prime order, and let `A_j` be the signed Weil block on `E_j`. Since the **whole** four-dimensional generalized inequality (19.13) is established and `M` is diagonal across these sources, minimizing over all old coefficients proves the genuine one-new-mode Schur residual

\[
\boxed{
D_j-B_j^*A_{j-1}^{-1}B_j
>\frac{35999}{10000}(1+p_j^2)G,
\qquad p_j=3,5,7.}
\tag{19.15}
\]

The first step adds `p=3` within the historically reviewed IV radius `11/10`; subsequent steps require the larger radii `17/10` and `2`. Thus (19.15) is an explicit finite sequence of mathematically valid original IV010-shaped **special-source** support extensions. It is not the actual all-old-core IV010/IV011 conclusion.

## 4. Arbitrary finite prime-set version (existence, not a uniform infinite family)

The same reasoning is not special to four primes. **For every fixed finite set `P` of distinct rational primes and every prescribed `C>0`**, one may choose a sufficiently small common `b>0`, and then `0<η\ll b`, such that:

1. all packets near `±log p` are mutually disjoint, with `a_p=log p+b+η`;
2. every off-diagonal original prime-power sample is zero, as no ratio `q/p` is an integer and no product `pq` is a prime power;
3. each complete original diagonal obeys `W(h_p,h_p)/||h_p||² ≥ −6+(3/5)log(1/[2(b+η)])`, which tends to `+∞` as `b→0`;
4. all normalized crossed Gamma entries are `O_P(b)` because mixed shifts remain bounded away from zero while `||h_p||_1=O_P(b)` and `G≥b−2η`.

Hence, for any prescribed `C>0`,

\[
\boxed{H_{W,P}> C\,M_P}
\]

for a suitably chosen `b(P,C),η(P,C)`.

This is a **finite-island theorem**, not a single fixed-width or fixed-radius positivity certificate uniform over infinitely many primes. The support radius must expand with `max P`, and the packet width may shrink rapidly as new prime-indexed interactions approach one another. No claim of positivity on arbitrary supported tests, or a basis/density argument capable of recovering such positivity, follows.

## 5. Exact F41/Phase IV receiving gate and governance

There are three different objects which cannot be identified without proof:

1. **This note's** `M_P`: the native physical **Gram of four newly constructed** smooth packet sources. This is diagonal by disjointness.
2. **This note's** `H_{W,P}`: the **complete original signed Weil** `4×4` matrix on those sources. The new nontrivial statement is `H_{W,P}>3.5999 M_P`.
3. **Historical F41's** `neutralNF606PhysicalH` (`242×242`) and `neutralNF606PhysicalV` (`242×4`): physical **L² integrals involving *different frozen* `nf545` and `NF587` generators**, not signed Weil matrix entries. Source blob above is directly pinned.

A transfer into the actual F41 bank would require explicit integrals `⟨nf545_i,h_p⟩`, source-defined selected-span containment, all relevant original signed mixed entries and certified residual/Schur data, plus exact preservation of the fixed historical numerical inputs. None has been computed or asserted here. The agreement in the number **4** with F41's mixed width `Fin 4` is purely coincidental and supplies **no** operator attachment.

In addition, the existing IV011 theorem needs a source-specific cross bound against **each** possible old Schwartz combination, not merely against the three chosen prime-packet directions. Arbitrary-support monotonicity and original Weil/RH positivity remain **unproved**.

### Certification matrix

| Obligation | Status |
|---|---|
| Live original IV010/IV011 source sign, physical correlation, prime cutoff | SOURCE INSPECTED |
| Historical original F41 physical H/V **type and fixed generators** | SOURCE INSPECTED; distinct from this H_W |
| Four genuinely nested supports, including original 11/10 radius | ANALYTICAL / EXACT AUDIT PASS |
| Both Mellin poles vanish for all four sources | EXACT |
| Original prime powers: diagonal `p²` only; every mixed prime side zero | EXACT / RATIONAL AUDIT PASS |
| Full mixed Gamma and diagonal Gamma constants | ANALYTICAL with exact rational thresholds |
| Signed original four-dimensional `H_W>35999/10000 M_P` | **ANALYTICAL PASS** |
| Three selected-source one-new-mode Schur residuals | **ANALYTICAL PASS** |
| Arbitrary finite distinct-prime set, suitably shrinking width | **FINITE-FAMILY EXISTENCE PROPOSITION** |
| Actual F41 242×242 H / 242×4 V attachment | NOT ESTABLISHED |
| Full infinite-dimensional old support Schur, operator closure, RH | OPEN |
| New Lean kernel / GitHub Actions certification | NOT PERFORMED |

**Independent reproduction:** `src/audits/cr_alp_019_prime_packet_family_audit.py`, standalone pure-stdlib `fractions.Fraction` (no float signs): verifies the literal support radii, the full prime-power inventory through `n=50`, strict positive/mixed logarithmic gaps, classic Gamma inequalities, and the four-dimensional row-sum/physical-metric Schur margin. Executed locally: **PASS**. The infinite-dimensional Gamma convolution estimates and the finite-set limiting argument are **human analytic proofs**; they are not Lean compiled by the rational checker.

**Next CR-ALP-020:** rather than claiming F41 attachment from an invented source basis, investigate whether an *actual historical `nf545OriginalSelectedRealGenerator`* has a provable mixed representer against one of these packet directions, including a projection into the genuine 242-dimensional physical selected span, with a strict residual Schur bound. If not tractable, document a source-specific **nonattachment** and keep the finite prime-packet theorem as a separate local result. **No writes to Weil-Lab, F41/F42, historical RH or provenance acceptance ledgers.**