# CR-ALP-003 — Smooth ramp admission, finite-table Gram carrier and genuine Lean interface

**Canonical-Return / CR-ALP research branch.**
**Opened:** 2026-11-22 PST (2026-11-23 UTC).
**Claim level:** DERIVED / ANALYTIC and finite regression PASS. The source-checked Lean declarations are imported **as source locators only**. **No new Lean kernel pass, no proof of (M), no new RH proof, no operator equivalence.**

## 0. Primary contracts re-opened for this stage

**Alpöge**, *A remark*, 20 Nov 2026, PDF printed p. 3, timeline items 4 and 5, https://alpo.ge/remark.pdf.

**Rubinstein et al.**, *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, working draft **v2, 22 Nov 2026**, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf:
- (I.2)–(I.4): the ramp \(S(x)\), centered error \(G(x)\), Mellin pole criterion.
- (II.13)–(II.18), p. 36–37: actual weighted prime-power table, centering \(M_{u,t,j}\) and the one-sided mean inequality (M).
- (II.21)–(II.24), p. 38: exact cell decomposition into remainders.
- (II.27)–(II.30), p. 39–40: character expansion and translation averaging.
- (II.33)–(II.36), p. 41: nondistinguished-cell error and reduction to \(G(x)\).
- Part III §III.8, §III.11: rank/trace input and theorem from twelve numbered hypotheses. This human-readable reconstruction **does not prove (M)**.

**Frozen Lean comparator:** public mirror at commit
7cbf7a0c32fca1badf05931879a37c3f8a2b03b4 of
https://github.com/CocoSgt/riemann-hypothesis-lean.
These are **source-verified anchors**, not independently compiled in CR. Relevant paths and declarations:

1. src/Formalism/SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean:
   definitions \`topExcess\`, \`GrowthInequalityHypotheses\`, \`SizeWindows\`, \`GrowthBoundAtPoint\`, \`GrowthInequalityDataExists\`.
2. src/Formalism/SizeWindowBrackets/ChainSlopeWindows/LipschitzWeightAndMinimum.lean:
   \`theorem proposition_at_point : GrowthBoundAtPoint\`; its proof explicitly invokes \`growthInequality\` and \`proposition_size_bound\`.
3. src/Formalism/SharpRestriction/GlobalTheoremInputs.lean:
   \`LinkG_at\`, \`LinkG\`, \`LinkS_at\`, \`LinkS\`, \`LinkDIV\`, \`GlobalTheoremInputs\`. \`LinkG_at\` includes the *local von Mangoldt class sum*, and \`LinkS_at\` includes Hilbert–Schmidt norm control. These are not definitions of the original signed Weil operator.
4. src/Formalism/RankGapBounds/MeanBoundSecondForm.lean:
   \`WeightedMeanData.Data.Main\`, \`DIV\`, and the slope/rank comparison data. This is **lattice/period/rank** geometry, not CR's physical \(L^2\) Gram pencil.

All cited Lean declarations were actually read at the frozen mirror ref, including the body of \`proposition_at_point\`. The analysis does **not** assert that these are identical to a Weil-positive spectral estimate.

## 1. A rigorous smooth-source recovery of CR-ALP-002's ramp

Let \(X>1\), \(L=\log X\), and
\[
f_X(u)=e^{u/2}\mathbf{1}_{[0,L]}(u),\qquad
k_X=f_X*\widetilde f_X
=\big(Xe^{-|y|/2}-e^{|y|/2}\big)_+.
\]
The previous phase computed the full signed Weil distribution \(\mathcal W(k_X)\), but \(f_X\) is not a smooth source. This phase closes the **original smooth-Schwartz source limit**, not the Connes pole-null trace admission.

Choose a fixed real, even, nonnegative \(\rho\in C_c^\infty(-1,1)\) with \(\int\rho=1\), and
\(\rho_\delta(u)=\delta^{-1}\rho(u/\delta)\) for \(0<\delta<1/4\). Put
\[
\boxed{f_{X,\delta}=\rho_\delta*f_X\in C_c^\infty(\mathbb R),\qquad
k_{X,\delta}=f_{X,\delta}*\widetilde f_{X,\delta}
=\sigma_\delta*k_X,\quad
\sigma_\delta=\rho_\delta*\widetilde\rho_\delta.}
\]
All kernels are real and autocorrelations. Since \(\sigma_\delta\) has mass 1 and support in \((-2\delta,2\delta)\), \(\mathrm{supp}\,k_{X,\delta}\subset[-L-2\delta,L+2\delta]\). The explicit ramp is globally Lipschitz with constant at most \(K_X=(X+1)/2\), so
\[
\|k_{X,\delta}-k_X\|_\infty\le 2K_X\delta.
\]
Both \(k_X\) and \(k_{X,\delta}\) have Lipschitz constant at most \(K_X\).

**Original quadratic form continuity.** Under the complete two-pole / signed prime / full Gamma convention of T1C and CR-ALP-002,
\[
\boxed{\mathcal W(k_{X,\delta})\longrightarrow\mathcal W(k_X)
\quad (\delta\downarrow0).}
\]
Proof components:

- **Prime sum:** all \(k_{X,\delta}\) have support in a common compact interval for \(\delta<1/4\), so only finitely many prime powers contribute; uniform convergence handles each term.
- **Poles:** dominated convergence under \(e^{\pm y/2}\) on the common compact support.
- **Full Gamma:** writing \(\Delta_\delta=k_{X,\delta}-k_X\), the only potential singularity is
\[
2\int_0^\infty[\Delta_\delta(0)e^{-t/2}-\Delta_\delta(t)]
\frac{e^{-t/2}}{1-e^{-2t}}dt.
\]
For \(0<t<2\delta\), split into
\(\Delta_\delta(0)(e^{-t/2}-1)+\Delta_\delta(0)-\Delta_\delta(t)\)
and use the uniform Lipschitz constants: the integrand is \(O_X(1)\).
For \(2\delta\le t\le1\), the numerator is \(O_X(\delta)\), and the weight is \(O(1/t)\), giving \(O_X(\delta\log(1/\delta))\).
For \(t\ge1\), exponential decay gives \(O_X(\delta)\). The Gamma diagonal term converges by the same uniform bound.
Thus
\[
\boxed{|\mathcal W(k_{X,\delta})-\mathcal W(k_X)|
=O_X(\delta[1+\lvert\log\delta\rvert]).}
\]
This is a **fixed-\(X\) continuity statement**: the implied constant is not uniform in \(X\), and therefore cannot be used alone for the \(\sqrt X\)-scale inequality.

The individual source can be translated by \(-L/2\) without changing its autocorrelation, and then has support in
\((-L/2-\delta,L/2+\delta)\). Therefore the source's required support radius grows at least as \(\frac12\log X+\delta\).

### 1a. Pole moments survive smoothing: a rigorous non-null attachment

Let \(m_\delta=\int\rho_\delta(v)e^{v/2}dv=\int\rho_\delta(v)e^{-v/2}dv>0\).
Then
\[
F_+(f_{X,\delta})=m_\delta(X-1),\qquad
F_-(f_{X,\delta})=m_\delta L,
\]
and
\[
\boxed{\mathcal P(k_{X,\delta})
=2m_\delta^2(X-1)L\longrightarrow2(X-1)L.}
\]
In particular **every** nonnegative symmetric mollifier preserves nonzero pole moments. This confirms the earlier pole-null obstruction and proves that the appropriate source route is *pole-preserving*, not pole-null.

**Standing:** complete original signed quadratic form on smooth compact source tests, by a genuine limiting argument. **Not** an operator-core theorem for the independently completed Suzuki or Connes carriers.

## 2. Exact finite-cell table as a family of physical Gram sources

Rubinstein (II.13)–(II.14) has for each cell a weight \(w_i>0\), a cutoff \(a>0\) and a positive-part contribution
\[
(a/r-w_i)_+.
\]

For \(a>w>0\), define the *cell source*
\[
f_{a,w}(u)=\sqrt w\,e^{u/2}\mathbf1_{[0,\log(a/w)]}(u),
\quad k_{a,w}=f_{a,w}*\widetilde f_{a,w}=w\,k_{a/w}.
\]
If \(a\le w\), take \(f_{a,w}=k_{a,w}=0\). Then, for **every** integer \(r\ge1\),
\[
\boxed{\frac{k_{a,w}(\log r)}{\sqrt r}=(a/r-w)_+.}
\tag{ALP3-cell}
\]
Its native physical mass and two pole moments are
\[
\boxed{\|f_{a,w}\|_2^2=a-w,\quad
F_+(f_{a,w})=\sqrt w\,(a/w-1),\quad
F_-(f_{a,w})=\sqrt w\log(a/w)}
\]
when \(a>w\). Correspondingly
\[
\boxed{\mathcal P(k_{a,w})=2(a-w)\log(a/w).}
\]
These formulas preserve the **nontrivial positive metric and both poles**. Each cell can be smoothed by the argument in §1.

For a finite set of cells \(I\), form the diagonal operator-valued autocorrelation
\[
K(y)=\operatorname{diag}_{i\in I}(k_{a,w_i}(y)).
\]
This is the autocorrelation of the diagonal \(L^2\)-valued source \(\operatorname{diag}(f_{a,w_i})\) and has positive-semidefinite Fourier transform **cell by cell**. Its native source metric is \(\operatorname{diag}((a-w_i)_+)\), which can be singular when \(a\le w_i\); it is **not** identified with the Hilbert-space metric of the theta-lift maps.

If \(i(r)\) labels the admissible cell of a prime power \(r=p^h\) and \(A(r)=E_{i(r)i(r)}\) is the corresponding diagonal selector, Rubinstein's raw table is the exact arithmetic sampling
\[
\boxed{
P_{u,t,j}
=\sum_{\substack{p\text{ odd}\\p\ne q}}\ \sum_{h\ge1}
(1-\chi_d(p^h))\,\Lambda(p^h)
\frac{\operatorname{tr}[A(p^h)K(\log(p^h))]}{\sqrt{p^h}},
}
\tag{ALP3-table}
\]
with \(a=a_t\) and the *actual* \(u,j\)-dependent \(w_i=w_{u,j}(i)\).

**Important:** this is an **exact lift of the arithmetic table into a matrix-valued *test-kernel***. It is not an equality with the original **scalar** zeta Weil form, because the place-dependent selector \(A(p^h)\) is essential. Omitting \(A\), character data, or individual cutoff weights changes the statistic.

## 3. Averaging is a finite character projection, not whole-operator positivity

For a finite abelian group \(H\) of order \(n\), \(T_2\subset H\), and \(\delta=|T_2|/n\), every column \(s\in H\) belongs to \(uT_2\) for exactly \(|T_2|\) translates \(u\). Hence
\[
\frac1n\sum_{u\in H}\mathbf1_{\{v(r)=1,\ t[r]\in uT_2\}}
=\frac\delta2(1+\chi_q(r))
\]
for prime powers where the quadratic sign is defined (Rubinstein (II.28), with his conventions for exceptional primes). Multiplying by \(1-\chi_d(r)\) gives, for prime powers coprime to \(2qd\),
\[
\boxed{\tfrac\delta2
[1+\chi_q(r)-\chi_d(r)-\chi_q(r)\chi_d(r)]}.
\tag{ALP3-characters}
\]
Thus translation averaging converts the labelled positive cell tests into the **signed combination of four character-twisted ramp sums** (Rubinstein (II.29)–(II.30)), after explicitly accounting for exceptional primes and primitive-character conventions.

The coefficient signs matter. A positive-semidefinite **diagonal test kernel** does not imply that the signed combination of completed Dirichlet \(L\)-function forms is positive. Each twisted explicit formula would require its own conductor, Gamma factors and possible poles. The averaging identity does not identify a physical source isometry or spectral operator conjugacy.

**Reproducibility:** the accompanying Python audit exercises the exact cell identity, finite-group translation counts, and four-character sign algebra; the local test run prints PASS. These are finite regressions, not certification of Rubinstein's (M) or the complete arithmetic construction.

## 4. Exact Lean-to-Weil boundary: intermediate theorem is now located

The source-checked intermediate Lean statement
\[
\text{\`proposition_at_point : GrowthBoundAtPoint\`}
\]
in \`LipschitzWeightAndMinimum.lean\` proves the implication (in the Lean source)
\[
\mathrm{GrowthInequalityDataExists}(C,\varepsilon,x)
\Longrightarrow G_z(x)\le x^{1/2+(C+3)\varepsilon}.
\]
The source definition of \(\mathrm{GrowthInequalityDataExists}\) **includes** \(\mathrm{topExcess}\le \mathrm{Err}\), character estimates, size windows, group data and per-cell weights. In particular, this is **not** a direct unconditional estimate of \(G_z\) from Weil positivity.

The actual difficult arithmetic input in the Lean development is upstream: one must manufacture data satisfying this package of hypotheses, especially its one-sided mean excess bound. The mirrored \`GlobalTheoremInputs.lean\` exposes further upstream named inputs \`LinkG_at\`, \`LinkS_at\`, \`LinkDIV_at\` for local prime valuations, Hilbert–Schmidt trace/rank control and lattice slope comparisons.

**No common object yet:** the Hilbert spaces and inner products used to prove the theta-lift rank/trace estimates are not the logarithmic-\(L^2\) source metric \(\|f_{a,w}\|_2^2=a-w\), and no specified source map intertwines these metrics and their form numerators. The common arithmetic statistic arises only **after** averaging, normalization, character decomposition and removal of other scales.

This is the specific point where a claimed "spectral projection of Alpöge's proof" would have to supply a new theorem, rather than simply rename an existing inequality.

## 5. Verdict and next test

**DERIVED:** source-ramp mollification and fixed-\(X\) whole-signed-Weil continuity; exact active/inactive cell autocorrelation dictionary; physical metric/pole formulas; elementary finite-group averaging and character signs.

**SOURCE PINNED:** the exact conditional Lean intermediate implication \`proposition_at_point\` and its named input definitions at frozen mirror commit. *Not* an independent local Lean compile.

**OBSTRUCTED:** direct two-pole-null CR-T1B trace admission for the nonzero-moment ramp; any claim that the averaged arithmetic estimate follows solely from positive definiteness of the source kernel.

**OPEN:** an \(X\)-uniform smoothing/control theorem with error smaller than \(X^{1/2+\varepsilon}\); completed operator-domain transport; independent proof of the arithmetic \(\mathrm{topExcess}\le\mathrm{Err}\) hypotheses; a physical-metric-preserving correspondence between the theta-lift trace/rank data and the Weil source Gram pencil.

**Next CR-ALP-004:** take Rubinstein (II.36), with the complete auxiliary character corrections, and write an *exact quantitative* chain from its mean inequality to the CR-ALP-002 renormalized Weil lower bound. Pin the upstream Lean theorem yielding \`GrowthInequalityDataExists\` and test whether any **noncircular** stronger source-level claim survives. If no physical carrier map exists, report **ANALYTIC REEXPRESSION, NO SPECTRAL TRANSFER**.

No changes to Weil-Lab, original RH/provenance or Hodge-RH custody.
