# CR-ALP-049 — Rational finite trial inverse and exact projected residual falsifier

DRAFT / NONPROMOTING. Conditional on the inherited original signed Weil continuum identities and pinned interval inputs. No Lean or RH result.

## Setup and trial inverse

Reuse CR-ALP-046 E85 interval bank SHA256 `d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22`, and CR-ALP-048 L=(72/25)I+Q85 P Q85 >=(51/200)I. H is the original low signed matrix in its nonorthogonal physical basis. Let B=Q85 T|E85. In each parity sector a midpoint inverse Cholesky provides a candidate V rounded to rational multiples of 10^-55. Direct 100-decimal outward interval congruence gives ||V^T H V-I|| <= max absolute row sum < 2*10^-54. Thus V is invertible; no numerical factorization is used for acceptance. Define B_V=BV. The high trial basis F uses e_j-(c_j/c_last)e_last for matching parity indices 86..205 (59 directions per parity). The matrix Y is obtained by rounding the midpoint solution of L_F Y=F*B_V to rational multiples of 10^-55; X=FY is exactly pole-null, with its coefficients enclosed by intervals. The original signed C=F*B_V and structured L_F are retained.

Let K=V^T H V-C^T Y-Y^T C+Y^T L_F Y. Exact outward midpoint-minus-row-radius Bareiss tests certify K>(1/4)I in 43 and 42 dimensions, respectively. The fitted in-band residual Gram row-sum upper bounds are 1.505e-97 and 1.459e-97. They refer to the chosen rational fit on 86..205, not to the full residual. Candidate generation uses mpmath and numpy; all entry tests use the inherited integer interval class. Saved rational coordinates and interval matrices are in result_1.json and result_2.json. The preliminary complete-tail envelope below is too conservative; it is not used to falsify any actual form.

## Full residual, including both pole bridges

In either parity let c_j=<e_j,b>, with b the genuine parity pole. For X in Q85, let u=R85(TV-(gamma I+P)X), where R85 is the ordinary sine complement and gamma=72/25. Then the residual R=B_V-LX equals Q85 u. For any real or complex column combination z, ||Rz||²=min_t sum_{j>85,matching parity}|(uz)_j-t c_j|². This identity includes every pole bridge. Split the tail into J=86..K and J'=K+1..infinity, with G_J=sum_J c², G_J'=sum_J' c², and means m_J=(sum_J c u)/G_J, m_J' analogously. Then

||Q85 u||² = ||P_F u_J||² + ||Q_K u_J'||² + G_J G_J'/(G_J+G_J') |m_J-m_J'|².

This is the exact orthogonal decomposition into the two independent pole-null blocks and the single bridge in this parity. The bridge follows by projecting onto the vector c_J/G_J - c_J'/G_J', whose squared norm is 1/G_J+1/G_J'. Both parity bridges are included by the global fit, but an isolated finite Galerkin computation cannot discard them.

## First all-tail bound, conservative only

For j>205 and source i<=85, v=nu_j>2u_i=2nu_i. The inherited sine formula h_ji=-(2/a)(v AI_i-u_i AI_j)/(v²-u_i²) permits pole subtraction -AI_i/(sqrt(a)b(-a))*c_j. Here |AI_j|<6, since v sum_m (1-(-1)^j exp(-2a beta_m))/(beta_m²+v²) <= 2/v+pi/2, and the total prime weight is <177/50. Expand the denominator around 1/v²; retain the signed first moment D=sum d_i nu_i, bounding the remaining terms by absolute cubic moments and weighted AI_i moments. The explicit resulting j^-2, j^-3, j^-4 envelope is implemented in audit.py. The boundary term has |<e_j,Tb>/c_j| <= 2 log j+60: the exact Laplace formula for J_beta b/c gives modulus <=4 v²/[beta(beta²+v²)] for v>=1, the decreasing resolvent sum gives g(v)-g0 <= 4 + (1/2)log(1+4v²) <= log v+5, |g0|<6, and the paired prime formula contributes <=4 sum omega when v>=2, |b'/b(-a)|<2 (even a larger bound 8 sum omega suffices). The term P X outside the band is bounded using ||P||<21/8; gamma X vanishes outside the band. All sums are dominated by reciprocal-power integrals and the decreasing log²(j)/j² integral beyond 205. Using the fitted in-band residual and columnwise Cauchy-Schwarz, the computed full residual Gram norm upper bounds are <9e22 (odd) and <3e19 (even). They are far too large; the very large energy-conditioned sine coefficients cause the absolute moments to lose cancellation. No useful full-tail acceptance follows.

## Exact finite-band rejection of the scalar residual correction for this trial

Extend actual signed coefficients and structured prime coefficients through index 410 using the same inherited interval formulas (24 Bernoulli terms, shift 512). The initial extension rejected an outward-rounded endpoint of atan's argument domain; use the exact atan addition formula with 1/4 and rerun. No change to source signs or pole constraints. The trial X remains supported through 205. With the original in-band pole fit, the matrix K-(200/51)G_fit,410 already has negative rational witnesses. This alone is not a falsifier of the actual projected residual comparison, since G_fit,410 is not necessarily a lower bound for R*R.

Instead use the exact finite pole-null orthogonal projection on F_{85,410}. Let U be the ordinary residual coefficient matrix through 410 (after any legal pole subtraction). Let c be the actual finite pole coefficient vector. The projected Gram is

G_proj=U^T U-(U^T c)(c^T U)/(c^T c).

Since F_{85,410} is a subspace of Q85, G_proj <= R*R. Thus a negative witness for K-(200/51)G_proj is also a negative witness for K-(200/51)R*R. It is not a negative original Weil witness and does not disprove H-B*L^-1 B>=0; the latter uses the true structured inverse residual R*L^-1 R, not the uniform scalar upper bound. The frozen rational coordinate vectors (denominator 10^12) verify the following strict upper bounds for their quadratic values: odd <-3, even <-2. The observed outward upper endpoints are -3.205209193370674 and -2.974708762284958. For the original chosen fit, corresponding upper endpoints are -4.441650539317954 and -4.111575746627835. The saved-matrix replay independently recomputes all 85 shifted Bareiss pivots and all four negative quadratic upper endpoints with integer arithmetic, without candidate optimization. It does not independently rebuild the continuum coefficients.

## Next gate

Do not merely refine the infinite remainder for this fixed trial and penalty: a finite projected lower bound already falsifies that sufficient test. Either extend the trial band substantially and retest the exact finite projected residual, or retain structured inverse information for the residual itself (e.g. a second Galerkin/residual correction, or a certified even-power Neumann remainder). A signed high-order moment expansion was considered but not executed once this finite falsifier was found. CR-ALP-048's 203-dimensional certificate remains valid conditional on its premises; the nonpositive-index bound 66 remains unchanged. Full collar, arbitrary support, Lean and RH remain open.

## Reproduction

Extract the archive with its directory layout, install numpy and mpmath if necessary, then run `python cr_alp_049/audit.py`, `python cr_alp_049/extend.py`, `python cr_alp_049/project.py`, and `python cr_alp_049/check_saved.py`. The last checker alone uses only the standard library and the saved result matrices. Each generator runs independently and pins the E85 input SHA256. The source interval formulas are inherited analytic premises, not kernel-certified. `extend_first_attempt.log` records the arctangent domain failure; successful logs and frozen matrices are included.
