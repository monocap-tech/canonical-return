# CR-INT-T1D-CC — matched joint-tail and directional error improvement

DRAFT / conditional on the inherited continuum/form identities and arithmetic-enclosure model. Executed on actual CA even 3-column and CB odd 12-column checkpoints. No new accepted dimension; conditional negative-index bound <=225 unchanged. No Lean, independent continuum review, full fixed-window positivity or RH.

## Joint tail
In the circular space, exact physical v and u are supported on I, exterior h on its complement. Thus full cross Grams v*h and u*h vanish. Let F denote the retained circular coefficients and T the omitted ones. With J=mI+D, ||D||<=r, the omitted sum -2 Re <Ju_T,v_T+h_T> is bounded above by -2m Re <u_T,v_T+h_T>+r(theta||v_T+h_T||²+theta^-1||u_T||²). The exact W_T Gram of v_T+h_T is G_v+G_h-(Z+Y)*(Z+Y) with directed errors and physical orthogonality, not the sum of two independent norm bounds. Omitted u*h midpoint is minus its retained unweighted cross; add it with sign and directed error to CB's base. Retain original B_u, J², finite exterior and signed energy terms. Exact Gram/action scale 10^12 and inherited rational interval definitions. Fixed theta,gamma: even (1/25,1/80); odd (1/50,1/175). Numerical parameter selection is heuristic, fixed verification uses exact rational arithmetic.

## Directional source error
Let Q=q*q/S², X=x*x/S^4 in the saved trial coordinates; let a=eY+sqrtup(tau), b=eU+sqrtup(tauU), e=representation error. For any positive lambdas summing to one, the squared action error is bounded by D2=a²Q/lambda1+b²X/lambda2+e²I/lambda3, by weighted Cauchy–Schwarz. Use lambda proportional to a||q||_upper, b||x||_upper, e. Mixed energy error is bounded by eK(eta Q+X/eta)+eD X, eta=||x||_upper/||q||_upper. Replace scalar allowances by these matrices and retain the same inverse Young gamma. Original physical Q Gram defect and conservative exact grid LDL retained.

## Executed fixed replay
| sector | inherited scalar/direct | separate matrix | joint | directional |
|---|---:|---:|---:|---:|
| even 3 | .2854024086714804 | .30246999999944585 | .310019999999432 | .3102199999994317 |
| odd 12 | CB .04506999999971855 | .04506999999971855 | .046069999999712306 | .04670999999970831 |

Odd CA scalar .089132 negative lower comparison and CB scalar .042980 negative lower comparison were earlier failures, not actual negative directions. Strict normalized rational margins: even >31021/100000, odd >467/10000. The exact shifted integer matrices at grid 10^12 are committed; positive principal minors verified locally. Fixed verifier uses no numerical eigensolver. The new directional bound does not add dimensions or establish new complement coercivity.

## Custody and reproduction
Use original six archives BP, BS, BT, BU, CA, CB. SHA-256 pinned in cc-inputs.sha256. Extract BP/BS/BT/BU flat into corresponding cr_int_t1d_<lowercase stage> sibling directories; CA/CB already contain named directories (CA includes odd). Place cc_fixed.py in their common parent. Run OPENBLAS_NUM_THREADS=4 python cc_fixed.py even and odd. Per-archive manifests verified before calculation. Outputs cc_results/even_terminal.json and odd_terminal.json reproduced fixed joint exact fractions and directed positivity for all four comparisons. Initial exploratory namespace collision (import overwrote src) repaired before executed tests. Original source data not regenerated; full chain custody remains open.

New terminal matrices and executed code are on this isolated CC branch to avoid racing the delegated archive upload on PR #2. Actual inputs still being transferred; no claim of clean GitHub-only source replay. Fourth even direction requires regenerating physical representation from BK even physical_dictionary.npz and exceptional.npz (plus inherited data). Those two are not in the local six uploaded archives. BK archive is listed as committed in chunks on the custody branch, but binary handoff to local execution remains unavailable. The fourth even trial has NOT executed. No new odd thirteenth trial. Next: recover BK physical dictionary and run four-column physical representation, directed circular action and matrix inverse. No new promotion until fixed four-column certificate passes.
