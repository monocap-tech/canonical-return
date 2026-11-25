# CR-ALP-063 — coalesced fourth signed-prime moment (candidate / nonpromoting)

Parent: 6673ef6605ecbf8529587b266650ce8ec5b40284. The exact frozen input archive remains SHA-256 f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2. No new witness selection, no finite high-sine cutoff. CR-ALP-061 completed with both four-scalar-moment comparisons intrinsically undecided.

## Physical computation

Write A=P_rat, Q=Q85, r_p for the projected physical proxy, and v0=A r_p-C for the exact midpoint-coordinate function constructed in 061. C is the finite sine and normalized pole-bridge correction, with exact binary128 coefficients; all exact basis functions are enclosed by the inherited MPFR interval arithmetic. The new evaluator reproduces the first two projections in the identical arithmetic order and the checker requires the displayed coordinate error eta to match 061 exactly.

Compute w=A v0=A² r_p-A C. On each half-domain cell, restrict the outer translations to arguments in [0,d] and the inner translations to final arguments in [0,d]. Exact rational knot tests ensure every active path belongs to a single source piece throughout the cell. Group paths by exact rational final displacement and source piece. Each group contributes one evaluation of r_p times the sum of products of its two positive prime weights. The finite correction contributes +sum of outer positive weights times C at the outer shifted arguments. This is precisely the original two negative signed translations, including the intermediate zero extension. No approximation is introduced by coalescing. The C++ independently checks that each stored path displacement interval overlaps the group displacement interval; equality and completeness are checked by the rational generator. The real and reflected domains have the same parity, so the half-domain integral is doubled.

The 061 complex-disk bound |v0|<1e42 applies for each active translated argument, by checking all final shifted arguments against the 058 disk geometry, h<=1/128, and the inherited bounds on the finite correction (which is entire within each active outer piece). Thus |w|<8e42<1e43. Products w² and w times normalized basis functions are bounded by 1e86; the exact degree-191 Taylor remainder bound with 96 Gauss nodes is <1e-12. Integrate ||w||², ordinary sine pairings and the pole pairing with outward intervals, then subtract the full finite orthogonal projection, exactly as in 059:
S_proxy=||Q w||²=||w||²-sum_i |<w,e_i>|²-|<w,b_perp>|²/||b_perp||².
This measures a complete physical norm, not a finite high-frequency sum.

## Transport and discrimination

The retained 061 displacement and projection arguments give ||Q P r-v0||<=delta=s eps+1e-13+eta, s=21/8, and ||(P-A)v0||<1e-12. Hence ||(Q P)^2 r-Q A v0||<=delta2=s delta+1e-12. If B is an exact rational upper bound for sqrt(S_proxy upper), then |S-S_proxy|<=2 B delta2+delta2², with S=||P_Q² r||². The checker computes B by integer square root and verifies overlap with the independent moment-support interval before intersecting the enclosures. It also checks that the reproduced second moment is consistent with 059 and that eta agrees with 061.

Use the four frozen rational linear weighted trials from the committed 061 feasibility audit. Evaluate the upper and lower inverse bounds over all 32 corners of the inherited G,M,N,T and new S boxes. Safe sufficient targets from 061 remain odd S<=1.147 or >=1.975, even S<=0.724 or >=1.467. The actual comparison will use the exact endpoint bounds, retaining prior bounds if tighter. Positive or negative signs concern this sufficient L comparison only, not a full original Weil sign. Full 85-column bank and fixed collar remain open.

## Custody

Workflow reconstructs the archive from Git, checks the predecessor data, repeats exact fourth-moment target generation and independently replays both 061 sectors before compilation. Both actual integration stdout records and exact endpoint JSON must be fetched back and committed before any result is accepted. No output left solely in CI artifacts. The new source is conditional on inherited classical analytical and forward-error arguments; no Lean theorem. Historical standalone archive reconciliation remains outstanding. Separate CR-ALP-060/062 matrix branch remains preserved.
