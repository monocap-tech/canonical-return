# CR-ALP-061 — pre-integration sharp decision regions (nonpromoting)

Inputs are the committed 058 and 059 exact rational endpoint JSON files, read at fb90b7ce3a55f931ede5984d8698b055db0bd979. No new physical moment is assumed. The native third-moment jobs were still running at inspection.

Let h=s²G-N, f=gamma G-M-ell u K, m=gamma G+M, d=gamma M+N, q=K-G²/m. For positive denominators and the nondegenerate four-moment formulas, sharp upper U3<K is equivalent to T>A, and sharp lower V3>K to T<B, where

A=gamma h+s²M-h²/f,
B=d²/m-gamma N+(M-dG/m)²/q.

This separates an intrinsic limitation of the four scalar moments from failure of a fixed cubic trial. Exact rational interval arithmetic on the inherited G,M,N,K boxes encloses A and B. Safe outward rounded display:

| Sector | A enclosure | B enclosure | Guaranteed undecided open interval |
|---|---|---|---|
| Odd | [0.051988,0.085653] | [-0.576190,-0.555945] | (-0.555945,0.051988) |
| Even | [0.252195,0.279856] | [-0.393059,-0.369359] | (-0.369359,0.252195) |

If the eventual entire T interval lies strictly in the last column and nondegeneracy is checked, both exact extremal inverse comparisons straddle K for every admissible tuple in the box. This is a sufficient preclassification only. For actual values use sharp_four_moments.py and also retain the previously certified fixed-parameter sufficient targets, which are tighter than these deliberately simple interval enclosures in some directions. Do not substitute these broad interval bounds for the fixed-parameter acceptance tests.

A separate necessary consistency check is obtained by taking Gram matrices of 1 and p against (s+p)dmu and (s-p)dmu, with p the centered spectral variable. This gives

(sM+N)²/(sG+M)-sN <= T <= sN-(sM-N)²/(sG-M).

The resulting conservative outer ranges are [-0.661833,0.462526] odd and [-0.632235,0.404366] even. These are not measurements. Any eventual interval wholly outside the corresponding range would require investigation of the numerical evidence or inherited premises.

The BigInt rational audit actually executed in the tool runtime on the fetched inputs and its exact result JSON are stored at src/audits/cr_alp_061/sharp_decision_regions.{js,json}. Replay with Node, passing the 058 verification_result.json and 059 second_result.json paths and optional output path. Both files fetched back and matched. No full moment bank, third moment, Lean theorem, or original Weil positivity established. Historical archive reconciliation remains open.
