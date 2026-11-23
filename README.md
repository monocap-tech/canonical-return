# Canonical Return

**Independent arithmetic–spectral comparison research archive.**

This public repository preserves the Canonical Return investigations from the signed Weil explicit formula through the comparisons with Suzuki, Morishita, Connes–Consani, and Deninger. Its purpose is to preserve exact source dictionaries, admissibility requirements, falsifiable transfer tests, and reproducible finite calculations.

**Status:** research notes / source-level analytical derivations / executable symbolic and rational-interval checks. **Not a proof of the Riemann hypothesis, not a proof of global Weil positivity, not a complete spectral equivalence, and not Lean-accepted formalization.** Some conclusions rest on published identities; other deductions have not received independent mathematical review.

## Research archive (complete Markdown snapshots)

The full text of the twelve research and migration Space pages has been archived under [`archive/space/`](archive/space/), in numerical order, preserving their canonical Markdown blocks. [`archive/SOURCES.md`](archive/SOURCES.md) records the exact page identifiers, original URLs and imported Git blob hashes. A historical migration-staging page is included for provenance; its statement that no remote existed reflects the *prior* staging date, not the current GitHub status.

| Sequence | Inquiry | Research page |
| --- | --- | --- |
| 00 | Charter and research jurisdiction | [Charter](archive/space/00-charter.md) |
| 01–03 | Suzuki: prime phase, scalar Weil form, support threshold | [SZ1](archive/space/01-suzuki-prime-kink.md) · [SZ2](archive/space/02-suzuki-full-weil.md) · [SZ3](archive/space/03-suzuki-threshold.md) |
| 04 | Morishita: orbitwise time reversal and global metric limitation | [MO1](archive/space/04-morishita.md) |
| 05 | Connes–Consani: source window and pole-null tests | [CC1](archive/space/05-connes-consani.md) |
| 06 | Deninger: Hodge-weight and affine time-reversal conditions | [DE1](archive/space/06-deninger.md) |
| 07 | Five-invariant comparison and formal compatibility model | [Integration](archive/space/07-integration.md) |
| 08–10 | T1A complex trace, T1B mixed Gram pencil, T1C interval enclosures | [T1A](archive/space/08-t1a.md) · [T1B](archive/space/09-t1b.md) · [T1C](archive/space/10-t1c.md) |

## New research (draft review branch)

[Draft PR #2](https://github.com/monocap-tech/canonical-return/pull/2) includes the source-level T1D sequence:

- [T1D-A](research/CR-INT-T1D-A.md): profile-uniform restricted positivity; first prime-29 collision.
- [T1D-B](research/CR-INT-T1D-B.md): prime-25 activation; wider two-cell positivity.
- [T1D-C](research/CR-INT-T1D-C.md): prime-8 two-entry perturbation and bounded width profile family.
- [T1D-D](research/CR-INT-T1D-D.md): prime-31/23/32 phase transport and exact complement Schur algebra.
- [T1D-E](research/CR-INT-T1D-E.md): original signed Weil mixed Riesz representers and rank-two projected penalty.
- [T1D-F](research/CR-INT-T1D-F.md): even/odd reduction and asymptotic Schur decoupling.
- [T1D-G](research/CR-INT-T1D-G.md): sharp-rate asymptotic residual law; exact prime-power arithmetic floor.
- [T1D-H](research/CR-INT-T1D-H.md): certified smooth-profile log-Fourier variance and nonzero finite-width mixed-coupling witnesses. [Checker](src/audits/canonical_return_t1d_h_mellin_probe.py).
- **[T1D-I](research/CR-INT-T1D-I.md): two-sided analytical finite-width projected residual bounds, \(3.2<\beta_+^2<6.2\), \(2.4<\beta_-^2<5.1\); resulting Schur penalty between \(320/1393\) and \(155/182\).** Coarse source-level enclosures, not high-precision quadrature, complement coercivity or RH. [Rational checker](src/audits/canonical_return_t1d_i_residual_bracket.py).

- [T1D-J](research/CR-INT-T1D-J.md): genuine two-dimensional translated/reflected smooth complement probes and a four-dimensional positive source compression. [Checker](src/audits/canonical_return_t1d_j_complement_probe.py).
- **[T1D-K](research/CR-INT-T1D-K.md): second reflected prime-4 probe family, a six-dimensional complete signed-Weil positive compression \(q>11/2\), and quantitative four-probe residual capture and Schur testing.** The outside-residual comparison additionally imports T1D-I; no whole-complement coercivity. [Checker](src/audits/canonical_return_t1d_k_extended_probes.py).
- **[T1D-L](research/CR-INT-T1D-L.md): explicit infinite-dimensional, very-high-frequency moment/Mellin-jet conditioned source with complete original Weil coercivity (>7) (at (B=2^{60})) and (>1800) (at (B=2^{2700})).** When importing the separate prior cross-block norm bounds, its authentic complement Schur condition holds on those deeply conditioned subspaces only; all unconditioned low-frequency and source-pole directions remain OPEN. [Checker](src/audits/canonical_return_t1d_l_frequency_tail.py).

- **[T1D-M](research/CR-INT-T1D-M.md): source-compatible log-prolate spectral conditioning and a reduced finite source Schur gate.** At \(B=2^{44}\), the original signed Weil form has a source-derived lower bound \(>15123/16000\) on an infinite-dimensional smooth pole-null subspace, with fewer than \(4B/3+951\) excluded prolate directions. Importing the separately derived T1D-I cross-block estimate gives a conditional Schur margin \(>136193/1456000\); **the finite low-mode block remains OPEN**. [Rational audit](src/audits/canonical_return_t1d_m_prolate_cutoff.py).

- **[T1D-N](research/CR-INT-T1D-N.md): new two-sided prolate rank estimate; at B=2^44, over 22 trillion genuine additional low-band constraints are necessary for T1D-M's uniform leakage condition, including a parity-sector lower bound.** [Exact-rational rank audit](src/audits/canonical_return_t1d_n_rank_barrier.py). This diagnoses the scale of the open gate; it does not solve it.


- **[T1D-O](research/CR-INT-T1D-O.md): exact physical-window prime-translation contraction** gives `||P_I|| < 55/4` (versus earlier `2 S_2 < 122/5`), lowering the admissible prolate cutoff from `B=2^44` to `B=2^29`. The complete signed Weil form is analytically `>164019/128000` on the corresponding infinite-dimensional conditioned smooth two-pole-null space. The further `>5005729/11648000` Schur margin **imports the separately unreviewed T1D-I residual bound**. The remaining roughly 683–716 million prolate low modes, full complement and pole directions are OPEN. [Rational audit](src/audits/canonical_return_t1d_o_window_prime_norm.py).

All new results are **analytic research assertions pending independent review**. Regression scripts check specified finite calculations and rational inequalities, not every functional-analytic step. Follow [Issue #1](https://github.com/monocap-tech/canonical-return/issues/1).

- [T1D-J](research/CR-INT-T1D-J.md): two authentic complement probes; four-dimensional restricted positivity.
- [T1D-K](research/CR-INT-T1D-K.md): six-dimensional restricted positivity and uncaptured mixed-coupling mass.
- [T1D-L](research/CR-INT-T1D-L.md): first source-level infinite-dimensional Mellin-jet-conditioned positivity.
- [T1D-M](research/CR-INT-T1D-M.md): logarithmic prolate conditioning and its finite low-mode Schur gate.
- [T1D-N](research/CR-INT-T1D-N.md): large necessary rank of the existing prolate leakage criterion.
- [T1D-O](research/CR-INT-T1D-O.md): support-compressed prime-shift norm bound <55/4; prolate cutoff 2^29.
- **[T1D-P](research/CR-INT-T1D-P.md): joint original prime-translation norm <93/10, verified by a rational three-step path sweep; prolate cutoff 2^22.** [Audit](src/audits/canonical_return_t1d_p_joint_prime_cubic.py). The whole low-mode Schur block remains unresolved.
- **[T1D-Q](research/CR-INT-T1D-Q.md): original joint prime operator norm <203/25=8.12 via exact-rational weighted Schur supersolution; source-conditioned cutoff B=2^21, still ~2.7 million excluded low prolate modes.** The conditioned cross-block margin additionally imports T1D-I; no full complement result. [Audit](src/audits/canonical_return_t1d_q_positive_supersolution.py).

- **[T1D-R](research/CR-INT-T1D-R.md):** exact original Gamma multiplier floor and sharper **negative** full pole-null semibound \(q>-2699\|f\|_2^2/200\), plus a **conditional** compact parity-separated Birman–Schwinger criterion for the authentic smooth complement. This is an exact spectral *test*, **not a certified positive low-mode eigenvalue**. [Audit](src/audits/canonical_return_t1d_r_compact_schur.py).

- **[T1D-S](research/CR-INT-T1D-S.md):** genuine parity Ritz **lower bounds** for the original compact Schur operators (even >7/25, odd >31/100), using 12 exact original smooth source cells; proves the old R.11 **generic tail estimator** requires an astronomical Dirichlet rank to become <1/2. **Not** a top-eigenvalue upper certificate or whole-complement positivity. [Audit](src/audits/canonical_return_t1d_s_ritz_tail.py).

- **[T1D-T](research/CR-INT-T1D-T.md):** exact **quadratic** low/high Schur-tail certificate exploiting that `6I` does not mix Gamma eigenmodes; conditional imported T1D-I bound reduces the half-unit *tail-estimator Gamma-energy threshold* from 5760 to 275. Two new even/odd broad-support smooth, doubly pole-null trial functions have bounded H¹ frequency scales and zero prime self-interaction. **No Ritz top-eigenvalue enclosure or full complement positivity.** [Audit](src/audits/canonical_return_t1d_t_tail_schur.py).

- **[T1D-U](research/CR-INT-T1D-U.md):** exact original seven-prime central-support image bound `||A P_J||<7/2`, and improved **true parity compact-spectrum lower witnesses** (even `>60/73`, odd `>300/377`) from sharp Gamma/Jensen estimates on the authentic broad smooth pole-null probes. [Exact rational/SymPy audit](src/audits/canonical_return_t1d_u_central_prime_gamma.py); [optional noncertifying FFT diagnostic](tools/diagnostic_t1d_u_fft.py). **No high/low Gamma projection estimate `eta_N`, eigenvalue upper bound below 1, or whole-complement positivity is established.**

- **[T1D-V](research/CR-INT-T1D-V.md):** two *independent genuine broad-support parity Ritz directions* in each sector, with exact rational H¹ mass/derivative Grams and Jensen Gamma bounds. They certify **at least two** positive compact Birman–Schwinger eigenvalues (even second `>600/781`, odd second `>3/4`), and sharper leading **lower witnesses** (even `>300/353`, odd `>600/731`). The prime self/cross Gram vanishes identically on the central four-space, but the true mollified Gamma and selected-source Schur Gram entries remain **uncomputed**. **No top eigenvalue upper enclosure, eta_N, whole-complement positivity, or RH.** [Audit](src/audits/canonical_return_t1d_v_central_ritz_multiplicity.py).

## Executable checks

The regression runner includes **26 executable audits** covering T1B, T1C, and T1D-A through W (including V-WIDE). Each has limited scope. T1D-I combines previously certified G/H rational intervals with a direct, finite-width Gamma-multiplier error estimate, off-support Fourier-tail control, exact prime-translation geometry, and the physical pole projection. This proves **coarse two-sided source residual bounds**. T1D-L/M additionally prove conditional-source coercivity on explicitly conditioned **infinite-dimensional** smooth subspaces, but **not** positivity on the unconditioned complement, the excluded finite block, or any completed operator.

Run with Python 3:

```bash
python3 -m pip install -r requirements.txt
bash tools/check_all.sh
```

These audits check **selected finite-dimensional assertions**, not trace-class domains, completed operator domains, or a full Hilbert-space spectrum. No automatic GitHub Actions workflow is installed.

## Research custody

See [`governance/CLAIMS.md`](governance/CLAIMS.md) for the claim levels, [`governance/WORKFLOW.md`](governance/WORKFLOW.md) for change discipline, and [`refs/primary-sources.md`](refs/primary-sources.md) for version-pinned literature. In particular, preserve both zeta poles, the signed prime and Gamma terms, complex pairing conventions, and the **native positive metric alongside the Weil-form Gram numerator**.

This repository is independent of `monocap-tech/weil`, `monocap-tech/weil-lab`, `monocap-tech/provenance`, and `monocap-tech/hodge-rh`. Archiving these investigations does not advance any of those repositories' accepted RH/Lean proof frontiers. No license is specified by this import.

- **[T1D-V-WIDE](research/CR-INT-T1D-V-WIDE.md)** (independent companion to the base V note): four authentic compact smooth pole-null test directions in **each** parity, retaining the original T1D-U direction. All internal prime-power Gram entries vanish exactly by support, while the original selected-space Schur rank-one coupling remains present. Exact rational unmollified H1 Gram and controlled mollification yield the stronger **true** compact-spectrum one-sided witnesses **even >8/9, odd >6/7**. The full mollified Gamma/cross-block matrix and spectral upper/tail are still open. [Checker](src/audits/canonical_return_t1d_v_wide_ritz.py). No Lean/global positivity/RH claim.
- **[T1D-W](research/CR-INT-T1D-W.md)**: actual broad-central **smooth** Gamma `2×2` matrices for `k=3,4` enclosed by a 2048-term exact-rational kernel series with a proven infinite tail and explicit mollification error; all 24 prime-power selected–central shifts leave **only 2 and 5**; the genuine parity Schur correction is uniformly bounded. Both Ritz eigenvalues of each authentic central two-mode pencil are **strictly below 1** (even `(120/151,37816/43155)`, odd `(10/13,105036/125125)`). These are *finite test-space* bounds, **not upper bounds on the full compact parity spectra**. [Audit](src/audits/canonical_return_t1d_w_full_central_gram.py).

- **[T1D-X](research/CR-INT-T1D-X.md):** original Weil central–satellite cross-coupling; complete prime 2/3/4 versus 7/8/9/11/13 census; positive true two-mode Weil determinants (>7/100 even, >1/5 odd), Gamma off-kernel rational bounds, and a nonzero second-stage Schur coupling conditional on T1D-I. **No full-complement positivity or RH.** [Exact-rational checker](src/audits/canonical_return_t1d_x_satellite_coupling.py).

- **[T1D-Y](research/CR-INT-T1D-Y.md):** direct original selected-to-satellite Weil bounds (all 24 prime powers checked; only 2,5,13,16,17,19 overlap), **without T1D-I's global residual estimate**. The actual rank-two selected Schur correction is bounded on both satellite directions and their central cross terms, establishing genuine Schur-corrected `2×2` determinants **>1/20 even, >1/5 odd**. The complete original Weil form is positive on the corresponding six-dimensional smooth pole-null space. **Not** positivity on the infinite-dimensional remainder, a top-spectrum upper bound, or RH. [Checker](src/audits/canonical_return_t1d_y_satellite_schur.py).
