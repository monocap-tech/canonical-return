# CR-ALP-123AD — Spectral receiving bound for the complete weighted response defect

123AC run 45933049682 has authenticated the frozen inputs and entered actual response-column evaluation. This is a staged post-receiver; no dependent run is authorized until the complete two-sector 123AC result has durable custody.

Let the weighted defect matrix be the map Z L sqrt(Delta) from the rational upper metric LDL decomposition of 123AC. Use the diagonal inverse majorant C^-1 <= (Q diag(d) Q)^-1 on the pole-null space, dropping the weighted pole subtraction only in the upper quadratic form. The resulting response Gram is bounded by the sum of the explicit prefix and infinite tail Gram. For both, replace trace by certified spectral norm bounds. No new physical integrations or response proposals needed.

## Explicit prefix

Use saved signed midpoint/radius coefficients and outward rational roots sqrt(Delta_l/d_j). Round each weighted coefficient into a rational box of grid 10^30. For midpoint matrix A and error E, ||E||_op <= ||E||_F. Thus the weighted prefix operator squared norm is bounded by (sqrt_up(lambda_max_upper(A* A))+sqrt_up(sum squared box radii))^2. Construct the midpoint Gram exactly with integer products. A symmetric grid rounding plus row-sum diagonal reserve upper-encloses that Gram. Certify its spectral bound by strict positive-shift tests using fraction-free elimination, not floating eigenvalues. Keep the better of this and 123AC's trace certificate.

## Infinite tail: nine-function common Gram

For each transformed response column, 123Y/AC preserved signed coefficients A,B,C_n and the residual norm. With J0=2047 and v_j=nu*j, the common leading basis consists of 1/j, log(j/J0)/j and seven cos(j theta_n)/j. In a given parity, the continuous main Gram of the first two functions equals [[1,1],[1,2]]/(2J0). Each cosine has main diagonal 1/(4J0). For the nonoscillatory 2x2 block, the difference between the discrete sum and half-integral is bounded entrywise by 3/J0^2: compare each sample to its preceding length-two cell, bound the error by total variation and the initial missing interval by sup/2. The functions 1/j^2, log(j/J0)/j^2, log^2(j/J0)/j^2 have total variation <= 1/J0^2 on [J0,infty); their sup <=1/J0^2. For cos^2 include half the same constant error and 1/(2 J0^2 |sin(2 theta)|). Abel bounds give nonoscillatory-cos cross <=1/(J0^2 |sin theta|) and cosine-cosine cross <=[1/|sin(theta_i-theta_j)|+1/|sin(theta_i+theta_j)|]/(2J0^2). All sine denominators are positively enclosed by original 046 arithmetic. The row sums of these entrywise error bounds give a diagonal reserve, hence an exact 9x9 common Gram upper B.

Let P be the rounded coefficient matrix with columns weighted by sqrt(Delta_l/dcap). The leading tail Gram is upper bounded by P* B P / nu_lower^2 plus the coefficient-box operator perturbation. The pole remainder and centered sine-series tail are charged separately using their weighted Frobenius upper bound. Thus tail spectral upper is the square of the sum of the leading spectral root, coefficient-box norm error, and remainder Frobenius root. Keep the smaller of this and the original tail trace bound.

Add the two spectral upper bounds to receive eta in 123AB. The finite b charge and pole trial U,E are unchanged. A positive final lower bound would certify the corresponding conditional pole sector only. If inconclusive, no negative witness follows.

## Verification scope

The script is staged with synthetic rank-one and diagonal spectral tests, exact matrix-order rounding and signed common-tail geometry. No actual 123AC output has been read, no numeric improvement asserted. The complete defect receipt must be independently authenticated by SHA256 before activating 123AD. No new runner.
