# CR-ALP-066 — full-frequency Gamma diagonal comparison on Q85

DRAFT / NONPROMOTING. Classical analytic derivation and exact rational arithmetic conditional on the inherited source identities CR-ALP-026/027/042/045/046. No Lean, no new witness acceptance, no complete collar, no RH. The exact Return-012 proof was not retrieved; its weights and prime estimate are NOT imported. Use the original signed source, a=11/10, L=2a, physical ordinary sine basis e_j, v_j=j*pi/(2a), seven negative prime shifts P, ||P||<=s=21/8, both original pole constraints. Q85 consists of ordinary high-sine coefficients j>=86 satisfying the corresponding parity pole constraint. The inherited layered proof gives Gamma>=gamma I, gamma=72/25, on this Q85.

## 1. Exact boundary correction

For each reflection parity let epsilon=(-1)^j, beta_m=2m+1/2, and f_beta(v)=v/(beta^2+v^2). The exact resolvent sine identity in CR-ALP-046 gives

  Gamma(q,q)=sum_j g(v_j)|c_j|^2 - B(q),
  B(q)=(2/a) sum_m (1-epsilon exp(-L beta_m)) |sum_j f_beta_m(v_j)c_j|^2.

Here g(v)=Re psi(1/4+iv/2)-log pi. To see the rank-one factor explicitly, K_beta e_j equals the bulk term 2 beta/(beta^2+v_j^2)e_j plus f_beta(v_j)/sqrt(a) times E_beta(t)=exp(-beta t)-epsilon exp(-beta(L-t)). Pairing E_beta with e_i gives 2(1-epsilon exp(-L beta)) f_beta(v_i)/sqrt(a). Subtracting K_beta in J_beta=2I/beta-K_beta produces the negative boundary correction. Opposite parities are orthogonal. The identity first holds for each finite beta partial sum by bounded K_beta and the sine expansion, hence for arbitrary L2 coefficients; the diagonal J_beta bulk terms are nonnegative. The full boundary series is bounded by the argument below. Monotone convergence of the nonnegative resolvent form and bulk series gives the identity on the logarithmic form domain, without asserting uniform finite sine convergence in that domain.

The exact pole moment says sum_j v_j/(v_j^2+1/4)c_j=0 (the omitted constant is nonzero for either parity). Thus in beta<32 terms we replace f_beta by f_beta-f_(1/2). Its absolute value is at most (32/120)^2 f_beta, since all v_j>120. This subtraction is legal for the whole constraint space, not just a fixed witness.

## 2. Infinite Schur bound

Use positive weights w_j=v_j^(-1/2) in weighted Cauchy-Schwarz:

  |sum_j f_j c_j|^2 <= (sum_i f_i w_i) sum_j f_j |c_j|^2/w_j.

For a centered small-beta row, use the absolute centered coefficients and the same bound, with two factors of (32/120)^2. In each parity the grid spacing h=pi/a>57/20 and v_j>=v0=120. For 16 small beta values (1/2,...,30+1/2), 1-epsilon exp(-L beta)<3/2. Since f_beta(v)<=1/v,

  sum_i f_i w_i <= v0^(-3/2)+2/(h sqrt(v0)),
  f_j/w_j <= 1/sqrt(v0).

Bound sqrt(v0)>109/10, and use 16*(3/2) for the summed weight. Multiplying by (2/a)*(32/120)^4 gives the exact small-beta contribution in the audit (<0.001319).

For beta>=32+1/2, the coefficient factor is <101/100. The nonnegative function sqrt(v)/(beta^2+v^2) is unimodal, with maximum <=beta^(-3/2). A grid sum of a nonnegative unimodal function is at most its integral divided by spacing plus its maximum: compare rectangles on both monotone sides, allowing one rectangle of maximal height at the peak. The integral on (0,infinity) is pi/sqrt(2) beta^(-1/2), by substitution and the elementary beta integral. Consequently

  sum_i f_i w_i <= [a/sqrt(2)+1/32] beta^(-1/2).

For fixed v, beta^(-1/2) v/(beta^2+v^2) is decreasing, so the step-two tail sum is at most one-half its integral over (0,infinity), namely pi/(2 sqrt(2) sqrt(v)). Thus after division by w_j the large-beta penalty is bounded by

  (2/a)*(101/100)*(a/sqrt(2)+1/32)*pi/(2 sqrt(2)).

All inequalities extend to the infinite sums by nonnegative summation and Cauchy-Schwarz. Rational replacements sqrt(2)>7/5, pi<22/7 give total penalty <1.685266<7/4. The exact BigInt arithmetic verifies Machin pi bounds, square-root inequalities, Taylor lower exponential bounds, both penalty budgets, and band evaluations. This is an analytical infinite-tail bound, not a cutoff on beta or j.

## 3. Source-preserving comparison

On Q85,

  Gamma >= diag(g(v_j)-7/4),
  W_Q >= Q85[ (1-theta) gamma I + theta diag(g(v_j)-7/4) + P ]Q85

for 0<=theta<=1, in quadratic-form order. For a rational computable diagonal use k_j=floor(j*(314159/100000)/(2a)), and D_j=log_lower(k_j)-37/20-2/k_j-7/4, using the original CR-ALP-026/045 staircase lower symbol estimate. For each band start one may extend its certified lower value across the band, because log(k)-2/k increases; use the exact rational log floor only at the lower endpoint. The same audit shows for j>=481 the diagonal lower exceeds gamma, j>=851 it exceeds gamma by >0.62085, j>=1251 by >1.00661, and j>=2049 by >1.50036. As j grows the actual D_j grows without bound. At j=86 the lower is only 1.18762, so the new diagonal alone does not dominate gamma on the entire Q85. No max-of-noncommuting-operators assumption is allowed.

With theta=1/10 the blended comparator remains globally coercive: subtract ||P||<=21/8, the lower scalar at j=86 exceeds 0.08576. The exact rational comparator can therefore be used in a finite-plus-infinite variational response solve without more scalar prime moments. For any trial x in its form domain with residual e=r-C_theta x in L2, the upper inverse quadratic is

  <r,C_theta^-1 r> = 2 Re <r,x> - <x,C_theta x> + <e,C_theta^-1 e>
                  <= 2 Re <r,x> - <x,C_theta x> + min_t sum_{j>85}|e_j-t k_j^pole|^2/d_j(theta),

where d_j(theta)=(1-theta)gamma+theta D_j-s>0 and k_j^pole is the appropriate true pole coefficient, not the integer floor k_j above. This follows from the diagonal lower bound and weighted Cauchy-Schwarz on the exact constraint space. Here C_theta denotes the compressed operator with rational diagonal D_j and P, not the full original Weil operator. Since W_Q>=C_theta, its inverse form is bounded by this response. The last equality/inequality requires x in the operator domain and full L2 residual, not only finite projected residual. The physical response can be assembled with a finite band, but all tail residual coefficients and pole constraint must be included.

## 4. Scope and custody

The exact audit was executed in a BigInt JS runtime and its source and output read back at 2edf4df2156a4221e2b63d5185d0ad375b04cf85. This verifies the rational constants; it is not independent proof checking of the classical analysis. The new lower comparison has a growing diagonal, but does not uniformly dominate the previous scalar-Gamma comparison L=gamma I+P_Q. No actual frozen even response bound has yet been evaluated against it. Next test: construct a source-pinned finite response for C_theta, use exact global residual coefficient enclosures and weighted pole-null inverse tail, compare against the frozen K. Avoid sixth-moment integration until this direct comparison is evaluated. Odd witness unchanged; prior 1048-dimensional positive subspace and conditional nonpositive-index cap 66 unchanged. Historical archive reconciliation open.
