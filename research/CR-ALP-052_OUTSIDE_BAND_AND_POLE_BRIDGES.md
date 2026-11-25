# CR-ALP-052 — Outside-band scalar residual falsifier and pole-bridge localization

DRAFT / NONPROMOTING. All continuum source identities are inherited classical premises. The previous 598-dimensional structured finite result remains valid; no new whole-collar positivity or RH result.

## Exact finite outside-band decomposition

Reuse CR-ALP-051's frozen degree-12 trial correction Z on F=span(e86,...,e600) intersect H0. In each parity write A=86..600, D=601..850 (matching parity), c_A,c_D genuine pole sine vectors, s_A=c_A^T c_A, s_D=c_D^T c_D. Let Y_A,Y_D be raw ordinary high-sine rows of TV-(gamma I+P)(X+Z), with gamma=72/25. Since Q85 removes the single relevant pole direction in each parity, projection onto a finite supported pole-null vector agrees with projection of these raw rows. The finite supported pole-null space through 850 decomposes orthogonally into F, the new pole-null band D, and one bridge per parity spanned by c_A/s_A-c_D/s_D. Its outside-F Gram is

G_out=Y_D^T(I-c_D c_D^T/s_D)Y_D + s_A s_D/(s_A+s_D) d^T d,
where d=c_A^T Y_A/s_A-c_D^T Y_D/s_D.

This is a lower bound on the full E*E, where E=(I-F)(B_V-L(X+Z)) lies in Q85. It does not assume that finite pole subtraction equals global pole subtraction.

## Executed exact arithmetic

Replayed CR-ALP-051's saved checker. Reconstructed Z using pinned 100-digit RF,D banks and 24-digit outward integer Chebyshev recurrence. Extended inherited actual signed coefficient generator through 850 at 100-digit scale, same shift 1024 and 24 Bernoulli terms. Reconstructed X from the pinned CR-ALP-049/050 rational rows, respecting exact final pole constraint. Corrected finite energy K'=K-2 RF^T Z+Z^T L_F Z uses exact symmetric RF^T Z. Both 43/42 shifted pivots pass (1/4 and 1/8). Raw new band rows include PZ before projection. The old band pole moment is computed using the exact sine prime matrix; no omitted Gamma term is substituted. Both frozen bank hashes and predecessor trial hashes checked. Midpoint eigenvectors choose rational candidates only, rounded to 10^-12; interval quadratic forms and independent standard-library replay certify acceptance.

Odd/even squared residual lower witnesses exceed 0.32643/0.32732 (targets 1/16 and 49/1600). Corrected K'-(200/51)G_out has normalized rational witness upper endpoints below -0.55336/-0.65808. Hence full scalar correction for this particular X+Z is false, even before adding the nonnegative in-band residual or further tail. On the norm witnesses, high-only contributions are about .0502643/.0927245 and bridge contributions .2761755/.2346032; exact bridge fractions exceed 84%/71%. This is not a negative original Weil direction and not a falsifier of K'-R'^*L^-1 R'. No full infinite upper bound attempted after exact rejection.

## Next: add both GLOBAL bridges rather than just extend F

Let c_in=S600 R85 b, c_inf=R600 b, with b the matching genuine pole function; s_in=||c_in||^2, s_inf=||c_inf||^2. Both are positive. The actual global bridge is proportional to c_in/s_in-c_inf/s_inf. It is pole-null and lies in Q85, orthogonal to F and Q600. The identity Q85 c_in=(s_inf c_in-s_in c_inf)/(s_in+s_inf) shows that the bridge belongs to E600, since Q85=Pi0 minus projection onto E85 and c_in is in S600. Thus Q85 H0=F orthogonal-sum span(two bridges) orthogonal-sum Q600 H0. Augmenting F by these two directions (not just finite-band bridge surrogates) gives exactly E600 minus E85. The norm constants are computable from pole norms minus finite sums. Each global bridge is a finite combination of genuine b and ordinary sines, so inherited sine-pole and pole-pole identities are relevant. Next calculate the original signed interaction with B_V and the structured L cross terms, then jointly re-solve the Galerkin system on F plus global bridges. The bridge-only residual fraction measured here is not a proof that such correction suffices. Exact scalar full inverse gate, full collar, Lean, arbitrary support and RH remain open. Nonpositive-index upper bound 66 unchanged.

## Replay

Archive extends the CR-ALP-051 directory layout. `python cr_alp_052/check_saved.py` verifies saved interval comparisons, 85 shifted pivots, both witness signs, exact fractions of bridge contribution and pinned bank hashes. `OPENBLAS_NUM_THREADS=1 python cr_alp_052/audit.py` regenerates interval matrices and candidates using NumPy/mpmath. The independent replay does not regenerate analytic coefficients or polynomial moments. `build.py` rebuilds audit.py and polynomial_core.py from the pinned predecessor scripts. Full finite result banks and logs are included. No accepting Lean theorem is claimed.
