# CR-INT-T1D-AS — Positive spectral surplus whitening by inverse compression

**Disposition:** derived analytical comparison and executed local directed-integer subspace certificates PASS in both parities. Conditional on AP's complete parity-circle comparison, AM's weighted arithmetic norm, AQ's certified deficit arrays and the explicitly stated AR binary64 integer-dot model. No independent review, Lean, remote CI, whole positivity, negative original witness, completed Suzuki–Connes transfer or RH.

## 1. Preserve the source and regain positive energy

Keep I=(-2,2), the two-pole-null smooth core, physical L² metric, AP's parity-specific zero-extension isometry E into the circumference-8 circle, all 24 original prime powers, hybrid weight 2/5, and AM norm 399/50. Let d_j be AP's hybrid lower symbol and t=1/2. AR and AQ exclude all j<1024 (512 constraints per parity) and use the certified upper normalized deficits b_j from AQ on the finite sets J. The inherited arrays contain 4628 even and 4691 odd indices; the prior scan covers every j>=1024, including its analytic infinite tail. This AS execution reuses the supplied, previously audited arrays; it does not rerun the 56-million-index deficit scan.

Let u_j=E*phi_j and T e_j=sqrt(b_j)u_j. For a finite set P of positive-surplus indices, independently certify s_k>=0 with d_k>=t+s_k. Since P and J are disjoint, and all remaining j>=1024 satisfy the inherited AQ comparisons, on the original low-coefficient-null source we have

    q0[f]>=t(<D_+ f,f>-||T*f||²),
    D_+=I+sum_(k in P) (s_k/t) u_k tensor u_k.

All unselected surplus is discarded, not approximated with an uncontrolled error. The constant I is Parseval on the original physical support. The pole constraints remain imposed; this comparison does not replace them.

## 2. Inverse-compression upper bound

On the parity circle let D=I+sum_(k in P)(s_k/t) phi_k tensor phi_k. Then D_+=E*DE>=I. For any physical g the variational formula gives

    <g,D_+^-1 g>=sup_f [2 Re <g,f>-<Ef,D Ef>]
                <=sup_z [2 Re <Eg,z>-<z,Dz>]
                =<Eg,D^-1 Eg>.

Since D^-1=I-sum_(k in P) [s_k/(t+s_k)] phi_k tensor phi_k, any rational alpha_k in [0,s_k/(t+s_k)] gives

    T*D_+^-1 T <= M-C C*,
    M=T*T, C_jk=sqrt(b_j alpha_k) <u_j,u_k>.

The upper matrix M-C C* is PSD, and its norm need not be small on all deficit coordinates. Suppose W is a subspace where its quadratic form is bounded by beta I, beta<1. Impose z=T*f in W (at most |J|-dim(W) constraints). Then by energy Cauchy–Schwarz,

    ||z||²=Re<f,Tz> <= <D_+ f,f>^(1/2) <z,T*D_+^-1 T z>^(1/2)
                      <=sqrt(beta)<D_+ f,f>^(1/2)||z||.

Hence q0[f]>=t(1-beta)<D_+ f,f>>=t(1-beta)||f||². The codimension argument yields the conditional negative-index bound |J|-dim(W)+512 per parity. The constraint space is not claimed to be a spectral subspace of the original operator. No previous mixed directions are subtracted.

## 3. Executed surplus certificate

`as_surplus.py` re-executes the unchanged AQ phase/log-table prelude and AM arithmetic norm proof. Propose the four neighbor offsets -6,-2,2,6 from every deficit frequency. Keep distinct indices >=1024 and <56,000,000 with certified d_lower>t, computed by exactly the same phase/log grid as AQ. Put s=d_lower-t, and alpha=floor(10^8*s/d_lower)/10^8. Verify parity, disjointness from the deficit set, and strict positive weights. Accepted columns: 12560 even and 12751 odd. All original arithmetic signs and the original Gamma lower comparator remain unchanged.

## 4. Fixed-point matrix and spectral check

`as_whiten.py` imports AR's interval matrix, exact-dot and row/column norm functions. Form M0 on the S=2^24 grid as in AR, with error eta_M. For each 512-column block of C, enclose square roots on the 10^8 grid and signed overlaps on the 10^9 grid, then form C0 on S grid by interval midpoint. Off-diagonal formula is exactly AQ/AR's (all indices positive, same parity, and disjoint). Accumulate exact absolute error row sums, maximum column sum; rho=ceil(sqrt(maxrow*maxcol))/S bounds ||C-C0||. Since max b<16 and alpha<=1, ||C||<=4 by Parseval and support compression. Therefore ||CC*-C0 C0*||<=rho(8+rho). The block products use AR's checked exact integer dot arithmetic: absolute dot product sums <2^53 by exact Cauchy bounds, assuming ordinary binary64 multiply/add or FMA without reduced precision. Accumulation across blocks is bounded by the sum of block maximum squared row norms, checked <2^63. Round the product down entrywise onto S grid, incurring at most |J|/S in operator norm. This gives an explicit total error for M-C C*:

    eta=eta_M+rho(8+rho)+|J|/S.

Use ordinary eigensolver only to propose W. Round eigenvectors/eigenvalues to S grid, verify AR's Gram/residual inequalities via exact integer dots and exact norm enclosures. Certify beta=(1+delta)(L+epsilon+eta)/(1-delta). No floating eigensolver result used as a sign certificate. Initial combined run passed even, completed odd whitened matrix, then timed out before odd spectral check; odd rerun separately with proposal threshold .993 instead of .995. The final script supports separate parity runs and uses those successful thresholds.

| parity | deficit count | surplus count | absorbed subspace | constraints incl. initial 512 | beta diagnostic |
|---|---:|---:|---:|---:|---:|
| even | 4628 | 12560 | 3765 | 1375 | .9972205985363597 |
| odd | 4691 | 12751 | 3823 | 1380 | .9971518098418853 |

Exact beta fractions:
- even 1410746592311722636797880585/1414678551949591842234826752
- odd 78996315630962480771141441145/79221954822995941847138828288

`as_check.py` reads the executed logs and verifies both exact fractions <499/500, verifies counts, and prints the combined conditional bound 2755. Thus q0[f]>=||f||²/1000 on the certified complement. This is 1282 fewer than AR's 4037 (31.76% reduction). Even and odd constraints sum to 2755; not a negative-eigenvalue census. Prior cutoff and mixed certificates not combined. No new negative original witness.

Reproduce (Python 3, NumPy, SciPy; standard library): optionally regenerate the inherited AQ arrays with `OPENBLAS_NUM_THREADS=1 python aq_scan.py`; then `python as_surplus.py`; `OPENBLAS_NUM_THREADS=4 python as_whiten.py 0 > whiten.log`; `OPENBLAS_NUM_THREADS=4 python as_whiten.py 1 > whiten_odd.log`; `python as_check.py`. Numerical proposals may vary and require rechecking; the log check records the particular executed certificates, not arbitrary reruns. The bundle includes the inherited arrays and dependencies. No remote CI or Lean.

## 5. Next gate — AT

The inverse-compression bound is still a relaxation: the variational supremum on the entire circle permits values outside I. Enforce the exterior zero condition in the dual variational problem, using a finite set of exterior test functions and an additional positive Schur subtraction. Alternatively introduce AP's exact positive Gamma image deficit and AM's arithmetic deficit into D_+. Remove the blanket 512 low-frequency constraints per parity only after evaluating their complete symbols and correct constant-mode normalization. Full fixed-window positivity, all-support positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.