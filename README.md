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

All new results are **analytic research assertions pending independent review**. Regression scripts check specified finite calculations and rational inequalities, not every functional-analytic step. Follow [Issue #1](https://github.com/monocap-tech/canonical-return/issues/1).

## Executable checks

The regression runner includes **eleven executable audits** covering T1B, T1C, and T1D-A through I. Each has limited scope. T1D-I combines previously certified G/H rational intervals with a direct, finite-width Gamma-multiplier error estimate, off-support Fourier-tail control, exact prime-translation geometry, and the physical pole projection. This proves **coarse two-sided source residual bounds**, not a completed operator or complementary positivity certificate.

Run with Python 3:

```bash
python3 -m pip install -r requirements.txt
bash tools/check_all.sh
```

These audits check **selected finite-dimensional assertions**, not trace-class domains, completed operator domains, or a full Hilbert-space spectrum. No automatic GitHub Actions workflow is installed.

## Research custody

See [`governance/CLAIMS.md`](governance/CLAIMS.md) for the claim levels, [`governance/WORKFLOW.md`](governance/WORKFLOW.md) for change discipline, and [`refs/primary-sources.md`](refs/primary-sources.md) for version-pinned literature. In particular, preserve both zeta poles, the signed prime and Gamma terms, complex pairing conventions, and the **native positive metric alongside the Weil-form Gram numerator**.

This repository is independent of `monocap-tech/weil`, `monocap-tech/weil-lab`, `monocap-tech/provenance`, and `monocap-tech/hodge-rh`. Archiving these investigations does not advance any of those repositories' accepted RH/Lean proof frontiers. No license is specified by this import.
