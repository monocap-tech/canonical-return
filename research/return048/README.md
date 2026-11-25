# RETURN048 — retained signed-form numerical audit

**Status: AUDIT-PASS for the archived interval-contraction replay; DERIVED numerical source claim pending independent review.** This is a self-contained custody snapshot of the RETURN047 parent and RETURN048 bundle. It does not promote an RH or global Weil-positivity claim.

## Scope and result

The producer `return048_retained.py` forms interval enclosures for the original retained signed energy and mass on two 43-dimensional real parity blocks selected from the RETURN047 receipt. It uses the pole, prime, and Gamma terms in its stated original signed form and constructs a rational congruence certificate. The saved receipt reports strict transformed diagonal-dominance margins and physical floors:

| Parity | Dimension | Transformed margin | Physical floor (approximate display) |
| --- | ---: | ---: | ---: |
| Even | 43 | 1 | 2.5066981220949833 × 10^-21 |
| Odd | 43 | 1 | 8.372330838477966 × 10^-19 |

These decimal floors are for orientation; the exact rational values and interval matrices are in `RETURN048_retained_receipt.json`. The producer's scope field explicitly reads `actual original retained signed form; no mixed complement gate`.

## Reproduction

From the repository root:

```bash
python3 research/return048/return048_replay.py
```

The replay needs only Python's standard library. It checks the SHA-256 of its RETURN047 parent, exact selected names and parity dimensions, interval symmetry and widths, integer congruence contraction, positive diagonal-dominance margins, mass row bound, Frobenius factor, and exact rational physical floors. It **does not** independently reevaluate the special functions or prove the producer's source identities.

The producer additionally requires `mpmath` and `numpy` and uses interval evaluation for its Gamma, prime, and pole formulas. Its original execution log is archived for custody; this import independently executed the replay, which printed PASS for both blocks. The analytic identity behind the producer and the interval error estimates require separate review.

## Custody

| File | SHA-256 |
| --- | --- |
| `../return047/RETURN047_payment_receipt.json` | `def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f` |
| `RETURN048_retained_receipt.json` | `469717b43ebcefb11cdcba90d1e73b39e076738ed042947468131c544d678833` |
| `return048_retained.py` | `eaaced4a5a9f5eadd5d428a057160aea717b0ed1daf77b8d3acd0a4f462f9834` |
| `return048_replay.py` | `a92f6192559c170d0414fe1e6728292822ff0e63c388c7af6614ff4e41860ffa` |
| `RETURN048_execution.log` | `d91411f5e6b330376afcf47c238aafc101d61971f7eee4489894b9bf39d2e505` |
| `RETURN048_replay.log` | `7e2bcfa7312d944925a52a701345fe75ad1d12832febd9b671004559b4241d51` |

Source bundle: `RETURN048_reproduction.zip`, SHA-256 `e77a457e22797234c2ff46a31c0d09d7f78ec2ede35ff9a7d7fd8b2faba7dfc3`. The six extracted members above are retained directly for review.

## Boundary

The two blocks establish a strict bound on the **selected retained subspace** as represented by this receipt. The mixed complement and its Schur transfer are unpaid. No full 242-source positivity, global Weil positivity, Lean certification, or RH conclusion follows from this snapshot. Do not change another project's acceptance ledger on the basis of this audit alone.
