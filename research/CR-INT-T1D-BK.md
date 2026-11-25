# CR-INT-T1D-BK — Explicit physical exceptional sector and direct Schur gate

**Status:** executed local directed-grid compact-basis and physical Gram audits PASS. Conditional negative-index bound <=240 (121 even +119 odd), now with explicit finite physical constraint maps. Complement relative coercivity margins >=17/10000 (even), >=53/10000 (odd). No original signed exceptional-sector Schur matrix has been evaluated. Conditional on the inherited original fixed-window (-2,2), smooth double Mellin-pole-null source identities, form closure/core transfer and ordinary binary64 exact-integer-dot arithmetic model. No new frequency scan, Lean, remote CI or independent mathematical review.

## 1. Compact exceptional basis

BJ saved pole vectors and residual corrections but not the accepted eigen-subspace basis. Reconstruct its enclosed central comparator Mhat and error eta from the actual BI matrices and BJ pole update at ridge 1/1024 even, 1/512 odd. Propose full eigenspaces via scipy.linalg.eigh, round the first k and remaining m-k columns to the 2^23 grid, call the resulting matrices Q,N (Q and N here denote the scaled rational matrices). Keep k=4960/5024, m=5081/5143. Independently verify delta=||Q*Q-I||, dN=||N*N-I||, cross=||Q*N|| and eps=||Mhat Q-Q Lambda|| by checked exact-grid products and row/column norm bounds. Both Gram errors <1. No reliance on eigensolver accuracy. Since Q and N have complementary dimensions, principal angles give ||(I-P_Q)P_{N-perp}||=||P_Q P_N||<=rho=cross/((1-delta)(1-dN)), using conservative inverse square-root bounds. For x perpendicular to N, write x=q+r with q=P_Q x and ||r||<=rho||x||. The restriction of Mhat to range Q has upper beta0=(1+delta)(||Lambda||+eps)/(1-delta); its mixed block norm <=eps/(1-delta). Thus the true constrained inverse comparator on N-perp has upper bound beta=max(beta0,0)+2rho eps/(1-delta)+||Mhat||rho^2+eta, where eta includes all BJ base, ridge, product, pole and constant errors. This avoids using a large global Lipschitz constant times the principal angle. Keep only N as compact rational coefficient basis; Q is reconstructed and checked when reproducing the audit.

| Parity | m | k | rank N | rho diagnostic | beta diagnostic |
|---|---:|---:|---:|---:|---:|
| Even |5081|4960|121|.0000343741|.9965717258689349|
| Odd |5143|5024|119|.0000336826|.9892962045311029|

Both exact rational beta <997/1000; additionally beta_even<1-2*(17/10000), beta_odd<1-2*(53/10000). Exact fractions in basis.log and exceptional.npz. No new exceptional dimension subtraction.

## 2. Actual physical constraint functionals

Let T be the synthesis of all certified deficit representers sqrt(b_j) u_j, with positive even cosine and odd sine circle modes u_j=cos(j pi x/8)/2 or sin(j pi x/8)/2 restricted to I, and even constant u_0=1/sqrt(8) appended last. b_j are the rational upper deficits in AQ/AW; dictionary stores exact numerator weights (denominator 10^8) and mode IDs in matrix order. Let g=cosh(x/2) or sinh(x/2), Pi=I-gg*/||g||^2. Define F=Pi T N and G=F*F. For f perpendicular to g, F*f=N*T*f. Hash-verify original AZ/BD physical deficit Gram M0, not the surplus-inverse comparator. On even augment M0 with true constant cross sqrt(b0/8) T*1 and diagonal b0/2; do not reuse BJ's residual constant cross. Compress M0 by N, then subtract the pole outer product of N*T*g, with inherited directed pole pairings. Products on 2^26 grid, all directed roundoff and input errors enclosed. Small symmetric eigenbasis is a proposal only; independently check Gram and invariance residual and use congruence to bound actual G from below and above.

| Parity | Gram error | lower diagnostic | upper diagnostic |
|---|---:|---:|---:|
| Even |.000062922354|1.771810146867|11.943137542963|
| Odd |.000096773860|1.879801988531|11.647887223897|

Exact checks certify (177/100)I<G_even<12I and (187/100)I<G_odd<12I. Thus F has full column rank; the physical complement has exactly codimension 121/119 in the corresponding pole-null Hilbert space. The constraint right inverse F G^-1 is bounded in L2 with norm squared <=1/lambda_min(G)<4/7 (since both lower bounds >7/4). This is a conditioning estimate for constraint coordinates, not the sign of q on their span.

On C=ker F* intersect g-perp, w=T*f is perpendicular to N. The constrained inverse bound gives ||w||^4=|<f,Tw>|^2 <= <Dphys f,f> <B_pole w,w> <= beta <Dphys f,f> ||w||^2. Here Dphys is the retained normalized positive energy (surplus plus original arithmetic slack), Dphys>=I. Hence q0[f]>=t(1-beta)<Dphys f,f> with t=1/2. Exact rational checks give parity margins above. The positive Gamma image energy remains available; it has not been used in this bound. The source-form completion and smooth core are inherited, not independently formalized.

## 3. Direct original physical Schur criterion (analytical, not a signed audit)

Let L be the original physical operator Gamma-A (compressed zero extension for the Fourier multiplier, signed 24 prime powers), and let C be the preceding physical complement in the pole-null Hilbert space. P_C=Pi-F G^-1 F*. Each F column is a finite trigonometric sum plus hyperbolic pole correction on I. Its zero extension has Fourier decay O(1/|xi|) by integration by parts, so the logarithmic Gamma multiplier times its transform is square integrable. Thus L F is an L2 matrix of representers. Inherited closed form and core hypotheses permit restriction to C; the restricted operator L_C associated with the closed form has inverse bounded by 1/mu_sigma. Any pole-null form-domain f decomposes uniquely f=F z+h with h in C, where z=G^-1 F*f. Since F lies in the form domain this decomposition preserves that domain.

Set A_F=F* L F, K=P_C L F. Completing the square in the closed form on C gives exact finite signed Schur matrix S=A_F-K* L_C^-1 K. Nonnegativity of S is equivalent to nonnegativity of the original fixed-window pole-null form (on completion, hence on core; converse uses inherited density). A negative direction of S yields a negative original-form direction Fz-L_C^-1 Kz in the completion and, by density, a smooth pole-null negative trial. No positivity follows merely from G>0.

For any proposed correction U with columns in C in the original operator domain, W=F+U has the same constraints as F. Put r=P_C L W. The Schur identity is invariant under this change of lift:

    W* L W - r* L_C^-1 r = S.

Therefore the explicit rigorous bracket is W*L W - mu_sigma^-1 r*r <= S <= W*L W. A sharper lower bound uses inverse of the compressed retained positive metric: r* L_C^-1 r <= mu_sigma^-1 r* (Dphys|_C)^-1 r. The scalar residual penalties can already be taken 10000/17 and 10000/53 instead of 2000/3. All matrices in this sign gate have size 121 or 119, but the operator response remains infinite-dimensional. Retain actual original Gamma boundary energy, prime terms and pole constraint when calculating W*L W and mixed residual. AF's near-null energy scales prohibit assuming a coarse absolute norm residual estimate will suffice. No original signed Schur sign, corrector solve or full positivity has been certified here.

## 4. Execution and custody

Two bk_basis.py runs and two physical_gram.py runs PASS; bk_check.py verifies 32 inherited input/manifest digest checks (with overlap), two additional original M hash checks in physical_gram.py, and exact recorded inequalities. No new carrier expansion. The selected basis is a numerical proposal verified by exact-grid products; no globally optimal basis claimed. The compact checkpoint stores actual integer N, exact physical dictionaries, enclosed small physical Gram, pole vectors/corrections and logs. It does not store full Q or upstream large inverse matrices. `prepare.py` generates BK scripts from BJ; existing saved scripts can be used without generation. Run from BI even/odd directories with `OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=1024 BJ_PARITY=0 BJ_OUTPUT=/absolute/path/cr_int_t1d_bk/even python /absolute/path/cr_int_t1d_bk/even/bk_basis.py` (odd use 512,1,odd). From the same BI directory run `OPENBLAS_NUM_THREADS=4 python /absolute/path/cr_int_t1d_bk/physical_gram.py 0` (or 1). Run `bk_check.py` from any directory. Upstream BI/BJ arrays, helpers, AZ/BD M grid and inverse checkpoints needed for full reproduction. No changes to upstream data. All rational diagnostics in executed logs; SHA manifest for new artifact files. The physical dictionary and N suffice to define the 240 exact functions independently of the eigensolver.

## Next BL

Construct original signed energy action on these explicit F columns or on better conditioned lifts with identical physical constraints, choose complement correctors U, and enclose the small signed Schur bracket with a frequency-dependent mixed residual bound. It may be preferable to choose smooth lift functions and exploit the exact Neumann Gamma boundary series, rather than directly evaluate high-frequency circle sums. A certificate must resolve sign at the actual near-null energy scale. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.