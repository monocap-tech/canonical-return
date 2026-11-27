# CR-ALP-123AE — Joint frequency Gram and exact pole-null subtraction

123AC run 45933049682 has authenticated its inputs and is evaluating its columns. This receiver is staged against the eventual complete two-sector receipt, not executed.

Let Z be the actual Q85-projected response-column defect, in the LDL-weighted frame (including sqrt(Delta)). For each column, sum_{j>85} c_j Z_j=0 exactly by construction. Set d_j=min(cap,gamma_lower_j)-21/8, dc=cap-21/8, with d_j=dc for j>2048. The full weighted pole-pair vector is therefore

    v = sum_{85<j<=2048} c_j Z_j (1/d_j-1/dc),

and the complete weighted pole mass is

    h = st/dc + sum_{85<j<=2048} c_j^2 (1/d_j-1/dc).

No infinite tail estimate for this inner product is necessary. The identities use real coefficients and fixed parity; the complex version uses adjoints. The source of st is the original full pole mass minus its low sine projection, not a truncated denominator.

The constrained inverse Gram is bounded above by

    G = sum_{j>85} Z_j^* Z_j / d_j - v^* v/h

by the weighted projection identity in 123F. In the receiver the rows denoted Z_j here already include the LDL weights; 123AC's stored prefix rows require multiplication by sqrt(Delta_l).

From 123AD, P is the exact rounded weighted prefix midpoint Gram, and Q is a positive upper Gram of the rounded leading nine-function tail. Their sum upper bounds the corresponding combined midpoint Gram. If a bounds its norm, and eps bounds the stacked operator perturbation (prefix box radius, tail coefficient radius, tail remainder), then actual unprojected Gram <= P+Q+(2 sqrt(a) eps+eps²)I. This follows by expanding the Gram of the actual stacked operator; different frequency blocks are orthogonal. It avoids adding separately maximized prefix and leading tail eigenvalues.

For midpoint m of the enclosed weighted pair v, error norm e and upper denominator hplus>0, the subtracted rank-one term satisfies

    v* v/h >= m* m/hplus - (2||m||e+e²)/hplus I.

Thus an exact spectral upper of P+Q-m*m/hplus, plus the two perturbation reserves, supplies another eta upper. Retain the minimum of this bound, unprojected joint Gram bound, separate spectral bound and original trace. All positive-shift arithmetic uses rational upward rounding. Final tiny pole energy roots use a 10^-110 grid, separate from coefficient root rounding.

The receiver includes synthetic checks for orthogonal rank-one Gram addition and final small-root scale. This does not imply numerical improvement on actual 123AC output. No new runner activated. Conditional complete pole signs remain open.
