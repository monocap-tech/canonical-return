# Canonical Return

For a provider-neutral restart, begin with [AGENTS.md](AGENTS.md). It points to the native governing ledgers and exact evidence checks; this README does not replace a live branch/run reconciliation.

**Independent arithmetic–spectral comparison research archive.**

This public repository preserves the Canonical Return investigations from the signed Weil explicit formula through the comparisons with Suzuki, Morishita, Connes–Consani, and Deninger. Its purpose is to preserve exact source dictionaries, admissibility requirements, falsifiable transfer tests, and reproducible finite calculations.

**Status:** research notes / source-level analytical derivations / executable symbolic and rational-interval checks. **Not a proof of the Riemann hypothesis, not a proof of global Weil positivity, not a complete spectral equivalence, and not Lean-accepted formalization.** Some conclusions rest on published identities; other deductions have not received independent mathematical review.

## Scientific disclosure program

**[Scientific Disclosure Charter — planned manuscript series and release gates](SCIENTIFIC_DISCLOSURE.md)**

A coordinated, cross-repository editorial initiative will document the mathematical research program, adaptive investigative governance, multi-agent architecture, verification, operational findings, originality assessment and reproducibility evidence. The charter establishes eight planned scientific deliverables in Markdown, LaTeX and PDF, with explicit claim boundaries and release review. It does not change the historical archive's standing or admit any new mathematical result. Manuscripts are not yet released.

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

## Executable checks

Two originally attached audit scripts are included verbatim:

- [`src/audits/canonical_return_t1b_exact_audit.py`](src/audits/canonical_return_t1b_exact_audit.py): exact SymPy arithmetic for the rank-two Mellin-pole nullspace, two-dimensional Hermitian Gram pencil and selected splitting.
- [`src/audits/canonical_return_t1c_rational_enclosures.py`](src/audits/canonical_return_t1c_rational_enclosures.py): exact-fraction outward analytical intervals for one explicit smooth test profile.

Run with Python 3:

```bash
python3 -m pip install -r requirements.txt
bash tools/check_all.sh
```

These audits check **selected finite-dimensional assertions**, not trace-class domains, completed operator domains, or a full Hilbert-space spectrum. No automatic GitHub Actions workflow is installed.

## Research custody

See [`governance/CLAIMS.md`](governance/CLAIMS.md) for the claim levels, [`governance/WORKFLOW.md`](governance/WORKFLOW.md) for change discipline, and [`refs/primary-sources.md`](refs/primary-sources.md) for version-pinned literature. In particular, preserve both zeta poles, the signed prime and Gamma terms, complex pairing conventions, and the **native positive metric alongside the Weil-form Gram numerator**.

This repository is independent of `monocap-tech/weil`, `monocap-tech/weil-lab`, `monocap-tech/provenance`, and `monocap-tech/hodge-rh`. Archiving these investigations does not advance any of those repositories' accepted RH/Lean proof frontiers. No license is specified by this import.
