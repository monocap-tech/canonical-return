# CR-INTEGRATION-001 — Five-invariant compatibility audit and a finite-place Hodge-weight model

**Independent Canonical Return / opened 2026-11-22 PST (2026-11-23 UTC).** Outcome: a typed cross-canonical compatibility matrix; exact analytic congruence and reflection lemmas; an explicit **formal model** combining finite-prime Mellin metric distortion, weight-one flow, an alternating pairing and quaternionic anti-linearity; a focused next theorem candidate. **Not Lean-certified, not a new arithmetic cohomology, not full spectral equivalence, not global Weil positivity or RH.** No source repository or stopped lineage was altered.

## 1. Sources actually compared

- [CR-SZ-002](https://chatgpt.com/space/page_8c63ddb9ebb481919002432bc72dab3b): Suzuki/Zeta23 signed pole–prime–Gamma equality and Fourier conventions on the smooth core, completed-domain identification open. [CR-SZ-003](https://chatgpt.com/space/page_f22d2c1da1f0819181937662f715a667): full two-bump phase/compression estimate, no complement bound. Primary [Suzuki 2606.09096v4](https://arxiv.org/html/2606.09096v4).

- [CR-MO-001](https://chatgpt.com/space/page_b6e89919038c8191bfbdd6d202544491): Morishita prime-circle anti-equivariance and orbitwise Haar unitary, not a global Hilbert or cohomology map. Primary [Morishita 2508.15971v5](https://arxiv.org/html/2508.15971v5), Theorem 3.6.

- [CR-CC-001](https://chatgpt.com/space/page_9bd218795c188191b89f97a8f592b814): exact \(u=e^x\) source/metric dictionary, pole-null four-cell phase calibration, CCM native metric distortion and conditional spectral Schur gate. Primary [Connes–Consani–Moscovici 2310.18423v2](https://arxiv.org/html/2310.18423), §4.6 equation (57), §4.7 Theorem 4.6; [Andrews–Connes–La Gatta 2611.13129v1](https://arxiv.org/abs/2611.13129) (Nov 9, 2026), conditional high-frequency/Schur route.

- [CR-DE-001](https://chatgpt.com/space/page_6d0eb966bbd48191b68024af7a4c4bef): Deninger’s *conjectural* weight-one positive cup/Hodge form and the affine reversal necessity. Primary [Deninger 2204.02714v2](https://arxiv.org/html/2204.02714v2), §2 properties 2.1–2.7, theorem 2.13.

- Original RH [PRED-SEAM-16](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_16_DENINGER_EQUALITY_20260914.json), [PRED-SEAM-17](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_17_SONIN_EQUALITY_20260914.json) and [RH-GLOBAL-STRUCTURE-8](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_GLOBAL_STRUCTURE_8_DORMANT_BRIDGE_QUANTIFIER_FILTER_20260921.md) retain STOP/carrier-transport obligations. The separate uploaded `Formal congruence transfer.txt` documents NF662's intended 15-generator real interval congruence; it is **a pending proof-program record**, not a certified theorem for this investigation.

## 2. Five-invariant comparison

| Source | Physical source/domain | Complete signed Weil pairing | Native metric / positivity | Arithmetic phases | Extension parameter |
|---|---|---|---|---|---|
| **Suzuki** | Compact smooth logarithmic tests; completed form domain distinct | **RECONSTRUCTED on common smooth core**: poles − prime + Gamma | Original physical \(L^2\) denominator / generalized \((G_a,K_a)\) | n9 two-bump phase retained, full mixed sign quantified | support radius \(a\); complement bound OPEN |
| **Morishita** | Prime-orbit circles and adelic/foliated phase spaces, **not** Weil Schwartz carrier | **NO** global Weil-form map | Orbitwise normalized Haar unitary; **NO** global Hodge metric | Fourier modes reverse; cross-orbit phase coherence OPEN | logarithmic dynamical flow reversed; not support enlargement |
| **Connes–Consani/CCM** | Multiplicative smooth window via \(u=e^x\), pole-null trace subdomain; separate Sonin space | **RECONSTRUCTED** scalar source; trace pairing only under its stated hypotheses | Multiplicative \(L^2(du/u)\) matches logarithmic mass; finite-place Sonin map is **nonunitary** without pullback correction | Four-cell pole-null phase fingerprint; real trace inputs first | window \(\lambda=e^a\) **and separately** finite place set \(S\) |
| **Deninger** | Conjectural arithmetic \(H^1\); no realized map from Schwartz | **NO** constructed source pairing equivalence | Positive Hodge inner product and quaternionic \(*^2=-I\) **conjectural** | anti-linear symmetry \(\rho\leftrightarrow\bar\rho\); not original prime mixed entries | weight-one cohomological flow \(\Theta\); not a proven support filtration |

**Status of arrows:** Zeta23↔Suzuki scalar smooth-core = PROVED analytically/RECOVERY (not novel; no new Lean). Zeta23↔CC multiplicative smooth scalar form = PROVED analytically/RECOVERY; target trace carrier = admissibility-constrained. Morishita↔adelic orbits = THEOREM/PROVED; Morishita↔CCM *global Hilbert* = UNKNOWN, not asserted. CC Sonin↔Deninger Hodge = NO CONSTRUCTED MAP; naive unshifted positive weight-one time reversal is algebraically BLOCKED. Suzuki↔Deninger cohomology = UNKNOWN. Original selected Picard/BAL carrier = NO TRANSFER.

## 3. Universal finite-dimensional rule: preserve the Gram **pencil**, not only a matrix

Let \(H=H^*\) be an exact Hermitian finite Weil-form matrix and \(G=G^*>0\) its *native positive physical metric*. For an **invertible source/basis map** \(T\), the transported objects are **both**
\[\boxed{H'=T^*HT,\qquad G'=T^*GT.}\]
Then \(\frac{z^*H'z}{z^*G'z}=\frac{(Tz)^*H(Tz)}{(Tz)^*G(Tz)}\); the generalized eigenvalues of \((H',G')\) are the same as those of \((H,G)\), and \(H\) and \(H'\) have the same inertia (Sylvester). Neither statement holds for an arbitrary noninvertible compression; a compression only tests a subset of directions.

**One exact counterexample to metric blindness:** \(H=\operatorname{diag}(1,4)\), \(G=I_2\), \(T=\operatorname{diag}(2,1)\) gives \(H'=\operatorname{diag}(4,4)\), \(G'=\operatorname{diag}(4,1)\). The *ordinary* spectrum of \(H'\) is \(\{4,4\}\), while its **correct generalized** spectrum remains \(\{1,4\}\). This explains why matching ordinary finite matrix entries or transporting only the numerator is insufficient for Suzuki/CCM comparison. The uploaded NF662 planned congruence program fits the algebraic pattern but cannot be promoted beyond its compiler/attachment status.

## 4. Exact finite-place reflection lemma — new-to-this-comparison *formal* model

Start with the **actual CCM Mellin multiplier** from equation (57), for a finite prime set \(S_f\):
\[\boxed{D_S(t)=\prod_{p\in S_f}\left(1-p^{-1/2-it}\right)},\qquad m_S(t)=|D_S(t)|^2.\]
Each factor has \(|1-p^{-1/2-it}|\in[1-p^{-1/2},1+p^{-1/2}]\); hence \(0<c_S^2\le m_S(t)\le C_S^2<\infty\) with \(c_S=\prod_p(1-p^{-1/2})\), \(C_S=\prod_p(1+p^{-1/2})\). Also \(D_S(-t)=\overline{D_S(t)}\) and **\(m_S(-t)=m_S(t)\)**. These are exact consequences of the CCM factor, not a claim about the full Sonin spaces.

On the **ambient** Mellin Hilbert space \(H=L^2(\mathbb R,m_S(t)dt)\oplus L^2(\mathbb R,m_S(t)dt)\), let \(h\) denote its positive Hermitian product, **linear in the first argument**. On the maximal multiplication domain set
\[\Theta(f,g)(t)=(\tfrac12+it)(f(t),g(t)).\]
Let \(Kf(t)=\overline{f(-t)}\) (antilinear involution), \(J(f,g)=(-Kg,Kf)\), and let \(R(f,g)(t)=(f(-t),g(-t))\) (linear reflection). Because \(m_S\) is even: \(K\) is antiunitary, \(J\) antiunitary with \(\boxed{J^2=-I}\), and \(\boxed{J\Theta=\Theta J}\) on the correct common invariant domain. Reflection \(R\) is unitary with \(\boxed{R\Theta R=I-\Theta}\). These properties express the **affine reversal** found in CR-DE-001.

Define the **explicit bilinear alternating pairing**
\[\boxed{\Omega((f,g),(u,v))=\int_\mathbb R m_S(t)\left[f(t)v(-t)-g(t)u(-t)\right]dt.}\]
Then, for all pairs in the common domain,
\[\Omega(x,y)=-\Omega(y,x),\quad \Omega(\Theta x,y)+\Omega(x,\Theta y)=\Omega(x,y),\quad \boxed{\Omega(x,Jy)=h(x,y)}.\]
The first relation uses the evenness of \(m_S\); the second uses \((\tfrac12+it)+(\tfrac12-it)=1\); the third is direct substitution. Thus the quadruple \((H,\Theta,\Omega,J)\) realizes the **formal algebraic properties** of a *positive weight-one, alternating, quaternionic* pairing with flow reversal. No zeta-zero spectral data are built into it: \(\operatorname{spec}\Theta=\frac12+i\mathbb R\) is continuous, and **this is not Deninger’s arithmetic cohomology nor a realization of Weil's form**. This is an elementary compatibility model, not an RH advance.

## 5. What the model actually separates

**Compatible algebra:** finite-place metric distortion, phase-sensitive complex Hilbert structure, the half-weight affine flow reversal and an antilinear \(J\) with \(J^2=-I\) coexist without contradiction. Thus **nonunitarity of the finite-place CCM map is not by itself an obstruction to a Deninger-type positive Hodge metric**; one must pull back the metric correctly.

**Still unconstructed:** a source-defined, injective/completed map from **actual complex pole-null Schwartz correlations** into this or another independently defined cohomological carrier that (i) preserves the *full signed* Weil form, (ii) respects prime-power phases and all Gamma/poles, (iii) intertwines appropriate generators/flow with the affine weight shift, (iv) realizes zeta zeros as the correct discrete spectrum, and (v) controls support and spectral tails. The toy shows **formal consistency only**. Deninger's 2023 real-coefficient no-go does not apply directly to this complex/quaternionic toy, but the toy does not meet the arithmetic functoriality requirements either.

**Parameter mismatch exposed exactly:** adjoining a new prime \(p\) to \(S_f\) changes \(m_S(t)\) by the nonconstant factor \(|1-p^{-1/2-it}|^2=1+p^{-1}-2p^{-1/2}\cos(t\log p)\). But for source tests supported in \((-a,a)\), the original Weil prime summand at \(p\) vanishes whenever \(p>e^{2a}\). Therefore **finite-place metric change** and **physical-support activation** are different operations; any commuting square must include the metric correction and the other terms. There is no unconditional equality of those two enlargements.

## 6. One focused *noncircular theorem candidate*: CR-INT-T1

**Target: a finite pole-null trace/Gram-pencil transport theorem for the original signed form.** Fix \(a>\tfrac32\log3\) and the [CR-CC-001](https://chatgpt.com/space/page_9bd218795c188191b89f97a8f592b814) explicit **four-cell two-Mellin-pole-null family**, initially its two independent real generators \(v_+,v_-\). Work only on their **complex span** after proving the target trace pairing's complex polarization and trace-class admissibility. Exhibit the **actual** Connes/CCM trace-test map and full trace expression, not a positive surrogate or a map defined by \(W\). Prove, entry by entry, on this finite domain:
\[\boxed{H^{\rm Weil}_{ij}=H^{\rm trace}_{ij},\quad \text{with both poles, full Gamma and all prime powers accounted for.}}\]
Separately prove the correct native metric matrix \(G^{\rm trace}\), its pullback by the actual source map \(T\), and the resulting **generalized-pencil congruence** if \(T\) is injective. Test the three unit phases using the already calculated exact distinct original prime contributions; any model that discards phase, uses the wrong pole-annihilation conditions or silently replaces the metric fails.

**Why noncircular:** The target trace map/operator is defined by the semilocal construction independently of the original Weil matrix and must satisfy independently established trace-class/domain hypotheses. No RH, Deninger Hodge existence, global positivity or spectral convergence assumption is permitted. The target is **finite and falsifiable**. A success would constitute a rigorous *restricted* cross-representation attachment; it would not directly prove Suzuki positivity, CCM's (H2')/(R1)/(R2), original STRICT/BAL, or RH.

**Precise first obligation:** inspect the 2026 trace theorem's hypothesis (real pole-null tests) and establish a valid **complex sesquilinear polarization** under actual semilocal operator trace-class control. If unavailable, document this as a definitive domain obstruction; do **not** assert full trace equality for \(\beta=i\). Once polarized, check the full 2×2 matrix, the metric correction and dependence on place set. This is a smaller, useful mathematical deliverable than asserting a global RH bridge.

## 7. Custody and verdict

**PROVED ANALYTICALLY HERE:** standard congruence invariance of Hermitian pencils and inertia; explicit finite-prime Mellin-factor parity; the doubled ambient Hilbert-space \((H,\Theta,\Omega,J,R)\) satisfies the *formal* Deninger-style weight-one/quaternionic axioms; place-set and support-radius changes are not identical. **CHECKED BY EXACT SYMBOLIC CALCULATION IN THE PRIOR CC NOTE:** the four-cell pole-null prime phase fingerprints (\(\beta=1,-1,i\)) are distinct; this is termwise prime information, not a full-form positivity claim.

**PROVED IN THE CANON, BUT NOT BY THIS NOTE:** Suzuki’s smooth-core operator identity, Morishita’s topological anti-equivariance, CCM’s finite-place Sonin isomorphism and Deninger’s restricted real-coefficient no-go. **CONDITIONAL/OPEN:** every global or completed-domain form intertwinement, positive arithmetic Hodge construction, spectral equivalence, semilocal complement estimate and original-RH endpoint transfer.

**Next:** CR-INT-T1A — check the precise Connes 2026 pole-null trace operator and its complex polarization/trace-class source assumptions; *only after that* seek a genuine two-dimensional mixed trace Gram equality. This independent Space remains a research notebook, not a certification or proof-status authority.
