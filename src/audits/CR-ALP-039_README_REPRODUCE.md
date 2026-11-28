# CR-ALP-039 reproduction and custody

**Research/nonpromoting only** — standalone exact arithmetic, not Lean or RH.

The source-pinned fixed physical collar is `a=11/10`, original source `Zeta23.EF.literatureRHS`, seven **negative** paired prime powers `2,3,4,5,7,8,9`, complete Gamma and two true pole moments. The reproduction holds the original 12-mode signed Hmid and physical M matrix bytes pinned to CR-ALP-034/038.

Requirements: Python 3.11+ standard library for both exact audits. Optional SciPy/NumPy is used **only** for the noncertifying signed 15,000-to-45,000-sine numerical pilot. No internet required.

From the extracted archive root:

```
python cr_alp_039/cr_alp_039_sharp_anchors_full_schur_audit.py
python cr_alp_039/cr_alp_039_independent_interval_ldl_verifier.py
sha256sum cr_alp_039/CR-ALP-039_full_twelve_shifted_H_interval_matrix.tsv
```

Expected exact shifted matrix SHA256: `a631f3cc7b99e47e4749a7931cd1307af336b5e20fe28b28946aa819ba56af77`.

The first program constructs the original signed 12-dimensional lowered Schur matrix, applies both analytically derived anchor upper bounds, and checks every shifted interval LDL pivot at `10^-25`. The second is an independently implemented directed integer interval LDL at `10^-20`, which rechecks the complete 12-mode matrix and exact anchor budgets from source constants.

**Analytical premises not formalized here:** digamma partial-fraction and series/integral variation inequalities, zero-extension physical translation kernel and sine coefficient formulas, original source prime cutoff, pole projector and source-domain identities, CR-ALP-027 uniform tail lower `1/4`, CR-ALP-034 signed Hmid PSD lower matrix, CR-ALP-038 interior source-specific residual bounds. The Python programs validate their exact rational consequences, not these continuum mathematical inputs.

**Not claimed:** positivity of the intervening 14,988-dimensional signed physical shell, entire 15,000-mode Schur PSD, original arbitrary-support positivity, native Lean compilation or RH.
