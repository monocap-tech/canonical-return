# CR-ALP-014 — Gamma-scale certification obstruction; sharper dimension-normalized cyclic convergence

**Canonical-Return / Alpöge projection, 2026-11-23 PST.**  
**Status:** SOURCE-PINNED / DERIVED finite-Hilbert inequalities / EXACT LOCAL AUDIT PASS. **No new Lean compilation, actual period-operator three-word instance, RH theorem, completed signed Weil metric/numerator isometry, or global positivity certificate.**

## 0. Frozen source basis

The primary external mathematical context is Alpöge, *A remark* (20 Nov 2026, https://alpo.ge/remark.pdf), and Rubinstein et al., **v2 working reconstruction dated 22 Nov 2026** (https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf). The exact formulas used here were directly reopened in the public Lean mirror, **CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4**:

- \`src/Formalism/LocalRingCardinality/DataParams₂.lean\`: \`Scale.L = x^eps\`, \`SizeData₂.CGam_pos\`.
- \`src/Formalism/SharpRestriction/FamCore.lean\`: \`CGam\` and \`CGam_eq\` fix a **positive** real constant along the family.
- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`: \`GamOf L=exp(L.CGam*m^5*f*L.sc.L)\`; \`LinkS_at\` bounds the full unprojected diagonal Gram defect and full distinct-word HS² overlaps **when its dimension-comparability antecedent holds**.
- \`src/Formalism/SharpRestriction/GlobalSchurBounds/Chosen.lean\`: \`GamL\` is the **same** exponential; \`ChProviso_at\` provides \(C_0^{-1}\le d_u\Gamma^2/D_*\le C_0\).
- \`src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean\`: \`RamRank_at\`, \`RamShareF\`, \`ramShareF_of\`; for nonzero rank and \(T_u>2,f\ge1\), \(q_u=d_u-\dim\mathrm{ramCu}_u\) is a **positive integer**, with \(q_u\le d_e d_u\), \(\sqrt{d_e}\le e^{-x^{\varepsilon/8}}\).
- \`src/Development/AbelianVarietyTheory/Definitions/AutomorphicMembership/EvenCeil.lean\`: the **specific** design choices \`mAt\`, \`nAt\`, \`fAt\`, not universal assumptions for every \`FamP\`, reinforce Gamma growth.

These are source statements and proof uses, **not new independent kernel validation**. No genuinely evaluated arithmetic 3-word \(KD\) instance or labelled nonzero cyclic phase-margin source certificate was extracted.

## 1. The new weak normalization improvement: \(O(\sigma^2)\)

Use three distinct words \(i,j,k\), \(d_u=\dim C_u\ge1\), full Gram \(K_{uv}=A_u^\dagger A_v\), ramified restriction \(\widehat K_{uv}=R_u K_{uv}R_v\), \(D_*>0,C_0\ge1,\Gamma>0\). Put
\[
\sigma=\exp(-x^{\varepsilon/8})\in(0,1],\qquad
\beta=\min(1,C_0\Gamma^2/D_*).
\]
Assume **simultaneously** the source hypotheses \`LinkS\`, \`ChosenProviso\`, \`RamShareF\` and the actual positive, proper \`RamRank_at\` regime. This means
\[
\|K_{uu}-I\|_{\rm HS}^2\le e_0d_u,\quad
q_u\le d_ed_u,\quad
\|K_{uv}\|_{\rm HS}^2\le(1+\eta_0)d_ud_v/D_*\ (u\ne v),
\]
\[
C_0^{-1}\le d_u\Gamma^2/D_*\le C_0,\qquad
\sqrt{e_0},\sqrt{d_e}\le\sigma,\quad \eta_0\le\sigma.
\]
Define the dimension-normalized cyclic trace
\[
\mathscr C_{ijk}(K)=
\frac{D_*}{d_id_jd_k}\operatorname{tr}(K_{ij}K_{jk}K_{ki}).
\]
CR-ALP-013 derived
\[
|\mathscr C(K)-\mathscr C(\widehat K)|
\le6(1+\sigma)\sigma\sqrt{\beta+\sigma\sqrt\beta}.
\tag{14.1}
\]
**New use of rank integrality:** \(1\le q_u\le d_ed_u\le\sigma^2d_u\) forces
\[
\boxed{d_u\ge\sigma^{-2}.}
\]
The **upper** side of \`ChProviso\` yields \(C_0\Gamma^2/D_*\le C_0^2/d_u\le C_0^2\sigma^2\). Hence \(\beta\le C_0^2\sigma^2\) and \(\sqrt\beta\le C_0\sigma\). Substitution sharpens (14.1) to
\[
\boxed{
|\mathscr C(K)-\mathscr C(\widehat K)|
\le6(1+\sigma)\sqrt{C_0^2+C_0}\,\sigma^2
\le12\sqrt{C_0^2+C_0}\,e^{-2x^{\varepsilon/8}}.}
\tag{14.2}
\]
This is an **absolute complex-number** error, not a relative/phase certificate.

## 2. The source's actual Gamma growth

Because \`Scale.L=x^eps\`, the fixed family constant \(c=\mathrm{Fm.CGam}>0\), and native ranks \(m\ge2,f\ge1\),
\[
\boxed{\Gamma(x)=\exp(c\,m(x)^5f(x)x^\varepsilon)
\ge\exp(32c\,x^\varepsilon).}
\tag{14.3}
\]
It follows for \(\varepsilon>0\) that
\[
\boxed{\Gamma(x)\sigma(x)\ge
\exp(32c\,x^\varepsilon-x^{\varepsilon/8})
\longrightarrow\infty.}
\tag{14.4}
\]
This uses only minimal source size conditions; the specific \`EvenCeil\` rank design grows faster but is not required.

The stronger **spectral-envelope** normalization is \(\mathscr B(K)=\sqrt{D_*}\,\mathscr C(K)\). The **specific** complement-based error upper expression of CR-ALP-013 was
\[
\mathcal U_B=6(1+\sigma)\sigma\sqrt{D_*}\sqrt{\beta+\sigma\sqrt\beta}.
\tag{14.5}
\]
Under \`ChProviso\` and \(d_u\ge1,C_0\ge1\), a case split \(\beta=1\) versus \(\beta=C_0\Gamma^2/D_*\) proves
\[
\sqrt{D_*\beta}\ge\Gamma/\sqrt{C_0},\qquad
\boxed{\mathcal U_B\ge6\sigma\Gamma/\sqrt{C_0}\to\infty.}
\tag{14.6}
\]
**What grows is an upper-bound *expression*, NOT the actual unknown error.** Even the improved weak estimate (14.2), when multiplied by \(\sqrt{D_*}\), produces an upper budget
\[
\mathcal U'_B=12\sqrt{C_0^2+C_0}\,\sigma^2\sqrt{D_*}
\ge12\sqrt{(C_0^2+C_0)/C_0}\,\sigma\Gamma\to\infty
\]
using \(D_*\ge d_u\Gamma^2/C_0\) and \(d_u\ge\sigma^{-2}\). Thus these **complement-derived certificates** cannot prove \(o(1)\) accuracy at the stronger natural scale.

An **independent, trivial** spectral envelope always holds:
\[
|\mathscr B(K)|,\ |\mathscr B(\widehat K)|
\le(1+\eta_0)^{3/2},\qquad
|\mathscr B(K)-\mathscr B(\widehat K)|
\le2(1+\sigma)^{3/2}.
\tag{14.7}
\]
So the best of the available bounds is uniformly \(O(1)\), **not divergent**. Neither source inequality proves a relative/phase error without a separately sourced **nonzero labelled cyclic margin**.

### Literal normalization mismatch with the fixed Weil triangle

\`ChProviso\` implies \(D_*\ge d_u\Gamma^2/C_0\ge\Gamma^2/C_0\to\infty\). By \`LinkS\`,
\[
|\mathscr C(K)|\le(1+\eta_0)^{3/2}/\sqrt{D_*}
\le\sqrt{C_0}(1+\sigma)^{3/2}/\Gamma\to0.
\]
CR-ALP-007's **fixed** three-cell signed Weil cyclic ratio \(\Xi\) is a **nonzero negative real constant**. Therefore the **literal equality** \(\mathscr C_{\mathrm{arith}}(x)=\Xi_{\mathrm{Weil,fixed}}\) cannot hold for all sufficiently large \(x\). This is a **normalization-specific** falsifier only; an \(x\)-dependent Weil source family, corrected scales and a genuine signed-form intertwiner are not excluded.

## 3. Exact finite diagnostic: error may shrink even when a certificate fails

For odd prime \(p\), choose \(f=2\), \(d=(p-1)^2\), \(r=(p-2)^2\), \(q=2p-3\), and **illustrative** \(\Gamma=p^2,D_*=d\Gamma^2,D_{\rm all}=3d,C_0=3\). Set full diagonal \(K_{uu}=I_d\); restricted three-word cross overlaps \(t/\Gamma\) with \(t=1/10\); complementary cross overlaps \((a,b,\pm ic)/\Gamma\), \(a=1/5,b=1/6,c=1/7\), using mutually orthogonal target sectors. The full positive source Gram is unitarily equivalent to \(F_{t/\Gamma}\oplus G_{\pm,a/\Gamma,b/\Gamma,c/\Gamma}\oplus I_{3(d-2)}\). Every source has the correct integer \(r,q\), the full finite-frame comparability ratio \(d\Gamma^2/D_*=1\), \(D_{\rm all}\le D_*\), and \(\|K_{uv}\|_{\rm HS}^2\le 1/(20\Gamma^2)<d^2/D_*\). Both orientations have identical complete spectra, HS² and restricted Bessel geometry.

Yet the full cyclic traces are
\[
\boxed{\mathcal T_\pm=\Gamma^{-3}(1/1000\pm i/210),
\qquad\widehat{\mathcal T}=\Gamma^{-3}/1000,}
\]
with opposite imaginary orientation for all \(p\), while the **actual stronger-scale absolute error** is
\[
\boxed{|\mathscr B(K)-\mathscr B(\widehat K)|=
1/(210d^{3/2})\to0.}
\]
This demonstrates that even **actual strong-scale absolute convergence** need not preserve phase when the reference cycle also vanishes. It also proves one cannot infer that actual errors diverge from the divergence of \(\mathcal U_B\).

**CRITICAL SCOPE:** These are exact abstract *finite* positive frames consistent with the displayed numerical rank/HS/comparability conditions; \(\Gamma=p^2\) is an illustration, **NOT** an instance of the source's exponential \`GamOf\` at a specified \(x\). No \`D7\`, \`BuRational\`, \`ChosenCode\`, period arithmetic, divisibility or original Weil form is instantiated.

## 4. Disposition / next

- **PASS / new derived theorem:** (14.2), weak normalized cycle error \(O_{C_0}(e^{-2x^{\varepsilon/8}})\), using *positive complementary rank integrality*.
- **PASS / source-specific scaling:** \(\Gamma\sigma\to\infty\), so the present stronger-scale completion error *certificates* are vacuous asymptotically; a different estimate or independent margin is necessary.
- **PASS / exact local regression:** standalone \`src/audits/cr_alp_014_scale_phase_audit.py\` with exact Fraction/SymPy checks, including prime rank-law samples \(p=3,5,11,101,1009\), native comparison numerics, Gram positivity, distinct source cyclic orientation, actual error, and scale-budget inequalities.
- **NOT ESTABLISHED:** any nonzero arithmetic labelled phase margin, numeric multiword theta-lift \(KD\), signed Weil \(H/M\) preservation, or global positivity theorem. **No new Lean compilation or CI claimed.**

**Next CR-ALP-015:** source-audit a genuinely attached labelled cyclic margin in arithmetic \`BesselId\`, \`Pu\`/period, lattice, or a higher-moment theorem, or explicitly close the *full operator correspondence* as **NOT ESTABLISHED / EVIDENCE INSUFFICIENT**. The conditional scalar arithmetic-to-smooth-original-Weil test-family reexpression from CR-ALP-004 survives separately. Keep Weil-Lab and other RH/historical ledgers unchanged.
