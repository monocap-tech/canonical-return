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

All new results are **analytic research assertions pending independent review**. Regression scripts check specified finite calculations and rational inequalities, not every functional-analytic step. Follow [Issue #1](https://github.com/monocap-tech/canonical-return/issues/1).

## Executable checks

The regression runner includes **sixteen executable audits** covering T1B, T1C, and T1D-A through N. Each has limited scope. T1D-I combines previously certified G/H rational intervals with a direct, finite-width Gamma-multiplier error estimate, off-support Fourier-tail control, exact prime-translation geometry, and the physical pole projection. This proves **coarse two-sided source residual bounds**. T1D-L/M additionally prove conditional-source coercivity on explicitly conditioned **infinite-dimensional** smooth subspaces, but **not** positivity on the unconditioned complement, the excluded finite block, or any completed operator.

Run with Python 3:

```bash
python3 -m pip install -r requirements.txt
bash tools/check_all.sh
```

These audits check **selected finite-dimensional assertions**, not trace-class domains, completed operator domains, or a full Hilbert-space spectrum. No automatic GitHub Actions workflow is installed.

## Research custody

See [`governance/CLAIMS.md`](governance/CLAIMS.md) for the claim levels, [`governance/WORKFLOW.md`](governance/WORKFLOW.md) for change discipline, and [`refs/primary-sources.md`](refs/primary-sources.md) for version-pinned literature. In particular, preserve both zeta poles, the signed prime and Gamma terms, complex pairing conventions, and the **native positive metric alongside the Weil-form Gram numerator**.

This repository is independent of `monocap-tech/weil`, `monocap-tech/weil-lab`, `monocap-tech/provenance`, and `monocap-tech/hodge-rh`. Archiving these investigations does not advance any of those repositories' accepted RH/Lean proof frontiers. No license is specified by this import.
