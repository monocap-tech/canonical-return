# CR-INT-T1D-AQ — Physical support-constrained spectral inertia

**Disposition:** derived analytical support-frame comparison; executed integer/Fraction audit PASS. Conditional on AP's original Gamma and all-prime representation and AM's arithmetic norm. No independent review, Lean, remote CI, full positivity or RH.

## 1. Source and positive baseline

Keep I=(-2,2), the original physical L2 metric, the smooth two-Mellin-pole-null core D^0, and the exact AP parity-dependent circumference-8 embeddings. In the even sector take j=0,2,4,... with normalized circle cosine modes (the constant mode has its own normalization); in the odd sector take j=1,3,5,... with antiperiodic sine modes. For positive j the normalized modes are cos(j*pi*x/8)/2 and sin(j*pi*x/8)/2 respectively. Write u_j for their restrictions to I and c_j=<f,u_j>. By Parseval for the zero-extended f, sum |c_j|^2=||f||^2 in each parity.

AP gives the original lower comparison q0[f]>=sum d_j |c_j|^2, where d_j=m_Gamma(j*pi/8)-(2/5)a(j*pi/8)-(3/5)(399/50), with the genuine symbol a containing all 24 original prime powers. The complete original Gamma image deficit and the residual arithmetic deficit remain nonnegative; they are not claimed to vanish. No physical support or Mellin condition is relaxed in the implication.

Take t=1/2. For j>=1024 let b_j be a certified nonnegative upper bound on (t-d_j)_+/t. Constrain c_j=0 for every j<1024 (512 per parity). The goal is to remove only some of the remaining positive-deficit indices, rather than all frequencies with d_j<0 or d_j<t.

## 2. Certified finite deficit census

`aq_scan.py` re-executes AM's exact arithmetic norm audit, and extends AP's integer-grid phase/log calculation to all 1024<=j<56,000,000. The AP rational Gamma floor m_Gamma>=log(j)+C at j>=1024 and the same global sum of coefficient upper bounds certify d_j>1/2 at j>=56,000,000. The diagnostic Gamma lower bound at this cutoff is 15.055326880026241; the sufficient full-symbol threshold is 15.0413167448. No finite or infinite frequency is omitted. The AP proof of the rational Taylor phase seeds and powers with error <1/4000 applies unchanged (both radix digits below 8192). The symbol upper slack is 1/150, and the log and Gamma terms are rounded downward on D=10^9 grid. All 24 weights and 48 signed shifts are retained. The only two arithmetic weights in the comparator are the fixed AP fractions 2/5 and 3/5.

For integer Gamma lower gam and symbol upper v on this grid, write the positive integer deficit r_j=max(0,2v+3(399/50)D+5tD-5gam). The upper bound b_j is ceil(r_j/25)/10^8, since t=1/2 and 5D/(2*10^8)=25. The scan gives 4628 even and 4691 odd positive deficit indices (excluding all j<1024). The NPZ arrays are derived data, not external premises; reproduce the scan before the support audit. Coverage is asserted sequentially.

## 3. Physical support Gram and weighted Schur test

For positive indices j,k of the same parity and j!=k, exact integration on I gives

    G_jk = (2/pi)[sin(pi(j-k)/4)/(j-k) + sigma sin(pi(j+k)/4)/(j+k)],

with sigma=+1 for cosine and -1 for sine. The diagonal is 1/2+sigma sin(pi*j/2)/(pi*j), so G_jj<=1/2+1/(3j). This is the Gram of the physically restricted circle modes, not of unconstrained circular Fourier coordinates. Since j±k are even, the sine in each off-diagonal term has absolute value exactly zero or one: one precisely when its argument index is 2 mod 4. Using 2/pi<637/1000, bound each nonzero term from above by 637/(1000*|j±k|). The audit uses the triangle inequality for the two terms; it does not use their possible signed cancellation. Each entry is rounded upward on C=10^6 grid, using integer ceilings. Independent Machin-series rational inequalities check the pi constants.

For a retained deficit set S define the finite-rank physical positive operator

    B_S = sum_(j in S) b_j u_j tensor u_j.

Its nonzero spectrum is that of the weighted Gram M=(sqrt(b_j b_k) G_jk)_(j,k in S), since B_S=TT* and M=T*T for the synthesis operator T. The positive Schur weights sqrt(b_j) yield

    ||B_S||=||M|| <= max_(j in S) sum_(k in S) b_k |G_jk|.

In particular, if that row bound is <999/1000, impose c_j=0 for all positive deficit indices not in S as well as j<1024. The exact AP comparator then gives

    q0[f] >= t||f||^2 - t sum_(j in S) b_j |c_j|^2
          >= t(1-999/1000)||f||^2 = ||f||^2/2000.

This argument does not require the removed functionals to be independent, and does not assert that the retained deficit frequencies are individually positive. It treats the retained deficits jointly in the actual physical metric. The finite constraints are compatible with the smooth pole-null core, on which codimension is no larger than their count. A negative-definite subspace of dimension exceeding that count would have a nonzero vector in the certified complement, a contradiction.

## 4. Executed certificate

`aq_support.py` initially removes indices whose single diagonal upper term exceeds the row budget. A greedy numerical score proposes additional deletions. All accepted row inequalities are checked independently on the final retained set with integer arithmetic, not inferred from the proposal or a floating-point eigensolver. Integer bounds prevent overflow in the final weighted products and sums. The final partition of positive deficit indices into retained and removed is verified exactly. Executed output:

| Parity | Positive deficits j>=1024 | Retained | Removed j>=1024 | Total constraints incl. 512 |
|---|---:|---:|---:|---:|
| Even | 4628 | 2753 | 1875 | 2387 |
| Odd | 4691 | 2784 | 1907 | 2419 |

Exact maximum directed row bounds:

- Even: 19969127399731/20000000000000 < 999/1000.
- Odd: 49796554296753/50000000000000 < 999/1000.

Thus the original smooth pole-null form is >=||f||^2/2000 on a subspace of codimension at most 4806, with conditional negative-index upper bounds 2387 even, 2419 odd, total 4806. This improves AP's conservative count of 6421 by 1615 (25.15%). It is a different exceptional system of coefficient constraints; no earlier mixed certificates are subtracted. It is not an actual count of negative eigenvalues, a global frequency cutoff, or an estimate of the sign on the remaining finite obstruction.

Reproduce from the bundle using Python 3, NumPy, and the standard library:

    python3 aq_scan.py
    python3 aq_support.py

Logs and derived NPZ arrays are included. Analytic AP/AM imports are explicitly conditional; the numerical audit does not prove the original explicit formula or resolvent comparison.

## 5. Next gate — AR

The next structural improvement is to keep signed physical Gram correlations rather than use absolute row sums, and optimize positive Schur weights jointly with the retained set. The norm target is a finite weighted support-concentration matrix with exact trigonometric overlap; its largest eigenvalue can be bounded by blockwise certified Schur complements or a rational positive supersolution for the entrywise majorant. Varying the baseline t requires a separately certified infinite tail and deficit scan. A separate route is to exploit AP's positive Gamma image correction and remaining arithmetic deficit in this physical metric, rather than reducing all positive energy to tI. None of these are yet certified. Complete fixed-window/global Weil positivity, Suzuki–Connes transfer, Phase IV and RH remain open.