# CR-MO-001 — Morishita: exact orbit transport, time reversal, and the missing Hilbert carrier

**Canonical Return / separate research lineage / opened 2026-11-22 PST (2026-11-23 UTC).** **Disposition:** CANONICAL RECOVERY + ANALYTIC ORBIT-LEVEL CONSEQUENCE; GLOBAL METRIC/WEIL TRANSFER OPEN. Not Lean-certified and not a new RH or general Weil positivity theorem. No modifications to Weil/Weil-Lab, original RH, or RH-HODGE repositories.

## 1. Source custody and exact scope

Primary: **Masanori Morishita**, *On a relation between Deninger’s foliated dynamical systems and Connes–Consani’s adelic spaces*, **arXiv:2508.15971v5 (2026-01-21)**, [source HTML](https://arxiv.org/html/2508.15971): §§1.1–1.2 (adelic \(\mathbb R_+\)-flow, prime orbits, Artin monodromy), §2.2 Theorem 2.2.9 (Deninger closed-orbit circles), §2.3 Theorem 2.3.2 (covering monodromy), §3 Lemmas 3.4 and 3.5, Theorem 3.6 (the comparison). This note uses **v5**, not an unstated earlier version.

Morishita's theorem constructs for finite abelian \(F/\mathbb Q\) a continuous map \(\Psi_F:\mathfrak X_F\to\mathscr X_F\) from Deninger-type foliated dynamical space to Connes–Consani adelic space. It is **Galois equivariant** and **\(\mathbb R_+\)-anti-equivariant**, and maps the indicated prime-orbit circles and unramified covering fibers onto the corresponding adelic objects. Its construction first restricts Witt-character data to roots of unity, passes to \(\check H\simeq\mathbb A^{\rm f}\), and finally maps \(((a_p),r)\mapsto((a_p),r^{-1})\) in §3; the final inversion causes the anti-flow.

In logarithmic time \(t=\log r\), if \(\Phi^D_t\) and \(\Phi^A_t\) denote the source/target flows, Morishita gives the **exact dynamical law** \(\boxed{\Psi_F\circ\Phi_t^D=\Phi_{-t}^A\circ\Psi_F}\) (where the maps and flows are defined). This law is topological: the paper does **not** assert a global unitary operator, Hodge metric, or original Weil pairing transport.

## 2. Positive result: the map is orbitwise an orientation-reversing circle isometry

For \(F=\mathbb Q\), Morishita §2.3 and Theorem 3.6 describe source prime circles \(\gamma_{p,\bar a}\simeq\mathbb R/(\log p)\mathbb Z\) that map **onto** the adelic circle \(C_p\simeq\mathbb R/(\log p)\mathbb Z\). Both are \(\mathbb R\)-orbits with the same **primitive** period \(L=\log p\). Choose an origin \(x_0\) in one source circle and origin \(y_0=\Psi_{\mathbb Q}(x_0)\) in its image. The anti-flow law gives the exact coordinate map \(\boxed{[t]_L\mapsto[-t]_L}\).

**Elementary proof of the orbitwise conclusion:** Every source point is \(\Phi_t^D x_0\); equivariance makes its image \(\Phi_{-t}^A y_0\). Since both stabilizers are \(L\mathbb Z\), equality of images is equivalent to equality modulo \(L\), so the restriction to this orbit is bijective and continuous with inverse represented by reversal. Thus it is a homeomorphism and an isometry for the **canonical orbit time metric**, with orientation reversed. This is stronger than merely stating that the *image* is the circle, but it follows immediately from Morishita's source theorem and the verified periods.

The same observation applies to **corresponding individual circles of length \(\log N\mathfrak p\)** in the unramified finite abelian coverings, whenever Theorem 3.6's image correspondence identifies those respective circles. This is only a local orbit statement: **no global injectivity, homeomorphism, or covering is asserted for \(\Psi_F\)**.

## 3. New-to-this-comparison functional-analytic calibration: canonical orbitwise Fourier transport

Equip each circle with its **intrinsic normalized Haar probability measure** \(dt/L\) (available without any global adelic or Deninger measure). Its \(L^2\)-spaces admit the exact map \(T_p:L^2(C_p,dt/L)\to L^2(\gamma_{p,\bar a},dt/L)\), \((T_p f)([t])=f([-t])\) after chosen origins. Then \(T_p\) is a **unitary linear isometry** (change of variable \(t\to-t\)); for arbitrary complex \(f,g\), \(\langle T_p f,T_p g\rangle=\langle f,g\rangle\). It is not the conjectural Hodge metric.

On Fourier characters \(e_m([t])=e^{2\pi i m t/L}\), \(\boxed{T_pe_m=e_{-m}}\). Thus the periodic-orbit derivative generator \(D=-i\partial_t\), on its natural periodic \(H^1\) domain, obeys \(\boxed{D_D T_p=-T_pD_A}\). Here \(D\) is **only the circle-time derivative**, not Deninger's conjectural global cohomological generator or Suzuki's full Weil operator.

**Without the paired origin choice**, the formula becomes \(T_p e_m=e^{i m\theta_p}e_{-m}\) for a unit-modulus origin-dependent factor. Consequently, the orbitwise Fourier **power/magnitude** and full complex inner products are invariant, whereas a *distinguished phase across different primes* needs coherent origins or further geometric data. This is a useful application of Weil-Lab's phase-fidelity discipline: an orbitwise unitary does **not** lose phase information, but it does **not** create a canonical cross-orbit phase comparison.

## 4. Why orbitwise unitary transport is NOT Weil-form or Hodge-metric transport

The first exact limitation is **type**. Morishita's map is between dynamical spaces; Weil-Lab's \(W_{\mathrm{original}}\) is a sesquilinear form on a specified Schwartz/source domain, retaining **both zeta poles, complete Gamma/digamma, and signed prime-power shifts**. A unitary equivalence between \(L^2\) of two **single closed orbits** has no automatic map from that Schwartz domain, no mixed prime-power matrix \(C(s,t;\pm\log n)\), no Gamma/pole operator, and no source-domain/completion theorem. Nor does a correspondence of periodic orbit lengths alone determine the coefficient \(\Lambda(p^k)/p^{k/2}=\log p\,p^{-k/2}\) in the Weil formula.

**Quantifier obstruction:** Morishita's theorem is about abelian number-field phase spaces and their orbit packets. The RH-facing test is a universal inequality over an infinite-dimensional space of complex Schwartz functions. Reproducing the set of prime-circle lengths does not imply the same quadratic form, its sign, or coercivity at every support radius.

**Metric transfer precondition (precisely typed):** If global measures \(\mu_D\) on \(\mathfrak X_F\) and \(\mu_A\) on \(\mathscr X_F\) were specified, an \(L^2\) pullback \(\Psi_F^*\) would be isometric on all target functions only if \((\Psi_F)_*\mu_D=\mu_A\), with the necessary measurability and domain hypotheses. A *surjective* global unitary would additionally require an essentially bijective measure-space correspondence, or a different proven range theorem. Neither global measure nor these conditions is established by Theorem 3.6. In noncommutative quotient settings the precise Hilbert construction itself must be separately defined.

**Weil/Hodge transport precondition:** Exhibit exact test-function/carrier maps \(S_D,S_A\) and independent forms \(Q_D,Q_A\) with \(Q_D(S_D f,S_D g)=Q_A(S_A f,S_A g)=W_{\rm original}(f,g)\) for admissible pairs; separately establish a positive Hodge pairing and compatibilities if invoking Deninger's speculative cohomological argument. A topological \(\Psi_F\) does not supply any of this.

## 5. Comparison against retained September RH analysis (no automatic re-entry)

The original [PRED-SEAM-16 audit](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_16_DENINGER_EQUALITY_20260914.json) already concluded: Morishita transports **phase-space/orbit** structures, not Deninger's positive Hodge metric or the selected Bombieri endpoint. The original [RH-GLOBAL-STRUCTURE-8 admission filter](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_GLOBAL_STRUCTURE_8_DORMANT_BRIDGE_QUANTIFIER_FILTER_20260921.md) requires an exact source substrate, carrier-preserving intertwiner, native metric, selected packet and support-uniform quantifiers. The new orbitwise unitary lemma operates on normalized Haar measure of **a single closed orbit only**. It does not meet those source, metric, and universal-quantifier gates.

There is no evidence here that STRICT marked-observability or BAL exact-gap rigidity is resolved or even directly bounded. Likewise, this comparison does not reactivate the stopped original RH-PRED or RH-HODGE lineage.

## 6. Classification and falsifiers

**RECOVERY (already canonical):** \(\Psi_F\) Galois equivariance, flow reversal, finite-prime orbit mapping, unramified monodromy compatibility.

**ANALYTIC CONSEQUENCE (elementary):** For each corresponding prime orbit with matched primitive period, \(\Psi\) restricts to an orientation-reversing circle homeomorphism; normalized orbitwise Haar pullback is unitary; Fourier mode numbers change sign, modulo origin-dependent phase. This is not presented as a new theorem by Morishita, or as research-level novelty.

**OPEN:** global measure-class preservation, full coherent phase across different closed orbits, actual Weil-Schwartz source maps, Gamma/pole/prime form equality, Deninger Hodge pairing transport, and any RH endpoint application.

**Falsifiers for alleged *global* transport:** a source/target domain mismatch; a pushforward measure differing from claimed target measure; a missing cross-prime sesquilinear entry; a sign/flow reversal error; an unaccounted \(k\log p\) repeat or von Mangoldt weight; no common operator/form domain; incorrect source selection or quantifier order. An arbitrary orbitwise positive pairing is not the original Weil pairing.

## 7. Next author and next gate

**CR-CC-001 — Connes/Consani.** Use Morishita's precise orbitwise reversal and Artin monodromy as **calibration data** while examining their semilocal Sonin/compression or adelic trace constructions. First ask for an actual continuous/test-function/representation map into the operator carrier with pairing and normalization—not merely matching prime-circle periods. Cross-check against PRED-SEAM-17's warning that trace nullspaces, Sonin space, and off-diagonal block kernels are different typed objects.

This author transition follows the user's proposed sequence and does not depend on establishing nonexistent global metric transfer from Morishita.
