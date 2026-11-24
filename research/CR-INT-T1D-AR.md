# CR-INT-T1D-AR — Signed physical support-Gram subspace certificate

**Disposition:** derived analytical finite-rank comparison; executed AQ deficit scan and AR fixed-point subspace certificate PASS. Conditional on AP's original parity-circle Gamma and arithmetic identities and AM's arithmetic norm. Exact integer matrix products are executed via ordinary binary64 dot products with checked exactness bounds (see §3); no trust in the floating-point eigensolver. No independent review, Lean, remote CI, whole positivity, Suzuki–Connes completion or RH.

## 1. Source and data

Retain the original physical interval I=(-2,2), the pole-null smooth source D^0, AP's parity embeddings, all 24 prime powers, hybrid coefficient 2/5, AM's norm 399/50, and AQ's baseline t=1/2. Let u_j be the normalized circular mode restricted to I. For j>0, even indices use cos(j*pi*x/8)/2 and odd indices use sin(j*pi*x/8)/2. AQ proves q0[f]>=t||f||²-t sum b_j |<f,u_j>|² on the subspace with j<1024 coefficients zero, where b_j are the nonnegative normalized upper deficits. The AQ scan was rerun to completion in this directory: 4628 even and 4691 odd positive deficits, all j>=1024 and the entire analytic tail covered. The arrays are derived data, not external premises.

The true signed weighted Gram is M= T*T, where T e_j=sqrt(b_j)u_j, so M_jk=sqrt(b_j b_k)G_jk. The exact Gram is AQ's cosine/sine overlap, with off-diagonal G_jk=(2/pi)[sin(pi(j-k)/4)/(j-k)+sigma sin(pi(j+k)/4)/(j+k)] and diagonal 1/2+sigma sin(pi*j/2)/(pi*j). Every sine used is exactly 0, 1 or -1 by integer residues. No absolute-value replacement of the off-diagonal entries is made.

## 2. Subspace certificate (not coordinate deletion)

Suppose W is a k-dimensional subspace of the m-dimensional deficit coordinate space, with M restricted in quadratic form to W bounded above by beta I for beta<1. Impose T*f in W, at most m-k linear physical constraints, in addition to the 512 initial low-frequency constraints in that parity. Writing z=T*f, we have

    ||z||² = <f,Tz> <= ||f|| sqrt(beta)||z||,

with absolute values or real parts as necessary. Thus ||T*f||²<=beta||f||² (including z=0) and q0[f]>=t(1-beta)||f||². The constraints are exact linear combinations of the physical coefficient functionals; W need not be an invariant subspace of M. This argument applies equally to complexifications of the real matrices. By the codimension argument on the original smooth pole-null source, n_-(q0)<=m-k+512 in that parity. No previously certified mixed dimensions are subtracted.

## 3. Rigorous fixed-point subspace check

`ar_signed.py` constructs symmetric rational M0=A/S with S=2^24, with certified ||M-M0||<=eta. The overlap entries use a Machin-series bracket for 2/pi, signed rational division, interval square roots of the exact b_j on a 10^8 grid, and outward integer multiplication and division. The error matrix is enclosed by the maximum absolute row sum of the per-entry errors. All interval products fit signed 64-bit integers by the explicit bounds in the code.

An ordinary numerical eigensolver proposes eigenvectors with eigenvalues <=0.997. Round those vectors to V=Q/S and eigenvalues to diagonal Lambda, also on the S grid. No proposed eigenvalue or vector is accepted without the following checks:

    delta >= ||V*V-I||,
    eps >= ||M0 V - V Lambda||,
    L >= ||Lambda||.

The products Q^T Q and A Q are integer matrices computed with ordinary IEEE binary64 dot products. For each product, exact integer squared row and column norm maxima a,b are computed after checking int64 sum-of-squares overflow bounds. The assertion a*b<2^106 proves by Cauchy–Schwarz that the sum of absolute integer products in every dot is <2^53. Therefore each individual product and every partial sum is exactly representable in binary64, irrespective of summation order in conventional dot-product GEMM. This uses the explicit arithmetic assumption of ordinary binary64 multiply/add (or fused multiply-add), not reduced-precision or approximate GEMM. The returned integer grids are checked for finite integral outputs; these checks alone would not suffice without the dot-product bound and arithmetic assumption. Residual and Gram differences are calculated on the integer S² grid. For any rectangular error matrix E, the exact bound sqrt(max absolute row sum * max absolute column sum), rounded upward with integer square root, bounds its operator norm. No eigenvalue convergence or approximate orthogonality is assumed.

Since delta<1, V has full column rank k. For any coordinate x,

    ||Vx||² >= (1-delta)||x||²,
    <M Vx,Vx> <= [L(1+delta)+(1+delta)eps+eta(1+delta)]||x||².

Here ||V||<=sqrt(1+delta)<=1+delta, and V* M0 V=(V*V)Lambda+V* residual. Consequently the rational quotient

    beta = [L(1+delta)+(1+delta)eps+eta(1+delta)]/(1-delta)

is a valid upper bound on M in quadratic form on W=Ran(V). This proves the desired subspace condition without a certified complete eigendecomposition, spectral gap, or positive pivot test of the entire signed matrix.

## 4. Executed results

| parity | m | k | constraints m-k+512 | eta (approx) | delta (approx) | eps (approx) | beta (upward rounded) |
|---|---:|---:|---:|---:|---:|---:|---:|
| even | 4628 | 3126 | 2014 | .000140071 | .000064217 | .000115190 | .996420 |
| odd | 4691 | 3180 | 2023 | .000279606 | .000065313 | .000119154 | .996903 |

Exact beta values from the successful execution:

- Even: 6578287973110657581429654903/6601922894669772650134372352.
- Odd: 19744399944025388440198539567/19805746988919178847710085120.

Both <999/1000, and a separate exact Fraction comparison of these logged results verifies both <997/1000. Thus q0[f]>=3||f||²/2000 on the certified complement, conditional on the imported original analytical identities and the stated dot-product arithmetic model. The negative-index upper bound is 4037 (2014 even, 2023 odd). AQ's 4806 is reduced by 769. This is an upper bound, not a negative eigenvalue census. The remaining constraints are linear combinations, not a global frequency cutoff. No negative original direction is certified.

Reproduction (Python 3, NumPy, SciPy, standard library):

    OPENBLAS_NUM_THREADS=1 python aq_scan.py
    OPENBLAS_NUM_THREADS=4 python ar_signed.py

The proposal may differ on another implementation. Every run must independently pass the integer-grid checks; the log is evidence only for the executed run. The AQ scan and signed check logs and reproduced arrays are supplied in the bundle.

## 5. Next gate — AS

Do not simply count additional signed Gram directions at a slightly different threshold. The reference-energy construction discards all positive surplus d_j-t and both AP positive deficits. Use the full positive diagonal background to whiten the physical support restriction, i.e. study T*D_+^{-1}T with D_+ the compressed positive energy plus a strictly positive floor. A finite certified resolvent lower/upper comparison must retain the physical support and pole constraints, and an upper enclosure must include the infinite tail. A smaller preliminary target is to eliminate the blanket j<1024 exclusion using certified low-frequency symbol bounds and correctly normalized j=0 mode. The present 4037 bound is not a positivity theorem. Complete fixed-window and global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.