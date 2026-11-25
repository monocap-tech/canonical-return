# CR-ALP-047 — Exact scalar Schur falsifier and actual signed finite band

DRAFT / NONPROMOTING. Conditional on the classical CR-ALP-046 analytic entry identities and pinned finite matrix. No Lean, global positivity, or RH claim.

## Source and finite band

The source remains the original signed Weil form on I=(-11/10,11/10), restricted to both Mellin pole moments. All seven original paired negative prime powers and the complete Gamma multiplier are retained. Let H be the pinned CR-ALP-046 E85 interval bank, SHA256 d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22. The basis order in each parity is interior w_j followed by anchor u_k. Let c_j=<e_j,b> with the matching genuine pole function. The exact high sine subspace F=span(e_86,...,e_105) intersect H0 is contained in Q85; each parity has 10 ordinary sines and one independent nonzero pole constraint. Thus F has dimension 18 and E85+F is a physically orthogonal direct sum of dimension 103. It is a proper subspace of E105, not the complete E105.

For high j and low interior i, W(e_j,w_i)=Hss(j,i)-(c_i/c_k)Hss(j,k); for the anchor W(e_j,u_k)=Hss(j,k)-(c_k/B)Hb(j). Here Hss, Hb, and B are precisely the analytic formulas of CR-ALP-046. For high j,l the same Hss formulas (including the complete trigamma diagonal) apply. The entire high band is represented by e_j-(c_j/c_l)e_l using the last high index l of each parity as pivot. Its actual signed block and all low-high cross entries are computed, without a uniform tail-energy replacement.

## Exact witnesses against scalar comparison

Write R(p,q)=|W(p,q)|^2/[W(p,p)||q||^2]. The scalar sufficient H85-(1/mu) B85 PSD implies R(p,q)<=mu for every p in E85 and q in Q85. This follows from Cauchy-Schwarz on the exact Q85-projected residual. In each parity use the saved midpoint only to search for a candidate and round its coordinates to integers over 10^55. The exact interval verification does not use the approximate Cholesky factor or eigensolver. The frozen two-mode witness uses q=e_89-(c_89/c_95)e_95 for odd indices, and the analogous 86,92 pair for even indices. Multi-mode witnesses have rational coefficients on all but the last high index; the last coefficient is defined exactly by the pole constraint. All candidates and certified endpoints are in result.json, integer endpoints at scale 10^100.

Odd pair R>53/100. Odd multi-mode R>105/100; even multi-mode R>27/100. Both multi-mode witnesses falsify the mu=1/4 sufficient criterion; the odd witness also falsifies the mu=1 scalar comparison. The actual W(q,q)/||q||² for the odd multi-mode witness lies near 3.5960706075, and |W(p,q)|²/[W(p,p)W(q,q)] near 0.2947080125. All four actual two-dimensional determinants are strictly positive, verified by intervals. These are not negative Weil witnesses. The scalar gate itself is rigorously false conditional on the source identities, rather than simply inconclusive because of arithmetic width.

## Complete actual finite mixed matrix

With the low H block, the exact high F block and both mixed transposes, the full 52-dimensional odd and 51-dimensional even interval matrices are symmetric by intersection of the two independent interval expressions. Compress outward to 65 decimals. Subtract the maximum row radius from the integer midpoint diagonal, then perform exact Bareiss leading determinant elimination with divisibility assertions. The respective row radii are 52 and 51 units of 10^-65. All 103 pivots pass. This proves positivity on E85 direct sum F, conditional on the analytic inputs. It does not include the two remaining finite E105 directions or any additional infinite complement.

## Arithmetic and reproduction

The coefficient generator adapts CR-ALP-046's interval primitives to scale 10^100, recurrence shift 256, and 24 Bernoulli terms; the explicit integral remainder bounds for psi and psi' are both <10^-92 (rational checks). Elementary functions use the inherited reduced-argument series, 400 terms, and checked error radii. Exponential sums use the same 256 terms and geometric tail bound. These are classical analytic assumptions, not formal Lean proofs. The saved E85 bank is checked by SHA-256 and compressed outward to the new scale. Run `python verify.py` in this directory with the pinned bank and result.json: exact frozen replay uses standard library only, reconstructs all new entries, compares every saved result field and passes all 103 pivots. `python audit.py` regenerates candidates using mpmath (85-digit midpoint search only); no acceptance depends on its approximate linear algebra. No independent derivation of the original analytic entry identities is claimed.

## Next obligation

Scalar coercivity mu=1/4 cannot close the E85/Q85 Schur gate, even with an exact infinite residual Gram. Keep actual tail energy in a weighted or block lower comparison, with rigorous control of interactions between the explicitly retained F block and the remaining infinite sector. A stronger diagonal comparison must be tested against the witnesses before investing in an entire coefficient sweep. The positive index reduction 66 from CR-ALP-045 is unchanged; full E85/Q85 positivity, complete collar, arbitrary support, Lean and RH remain open. The full infinite residual matrix was not constructed because the exact finite witnesses already falsify the proposed scalar gate.
