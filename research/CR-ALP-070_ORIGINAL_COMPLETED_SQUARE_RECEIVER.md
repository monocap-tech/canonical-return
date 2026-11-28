# CR-ALP-070 — exact original completed-square receiver

DRAFT / NONPROMOTING. The CR-ALP-069 frozen even residual comparison passed, but is not itself the original corrected Schur inequality. Original correction X denotes the *CR-ALP-056* complete frozen correction contracted against the normalized even K1250 witness. Response x denotes the distinct, recovered CR-ALP-067 finite response trial. Do not interchange them. Retain the original signed source, Q=Q85, L=gamma I+P_Q, gamma=72/25, ell=gamma-21/8=51/200. C=Q[D_new+P]Q=L+QAQ with A=D_new-gamma nonnegative diagonal from CR-ALP-068; A denotes the ordinary high-sine multiplier (and QAQ its compression). C_old=Q[D_old+P]Q, the weighted 067 comparator. All diagonals rounded down to 1e-9 per band, with final band constant. C>=Q diag(d_j) Q, d_j=D_new,j-21/8>0.

## 1. Receiver identity

Use R and H for the original low-tail coupling and low energy evaluated on the same witness. K_L=H-2 Re<R,X>+<X,LX>, r_L=R-LX. Complete the square against C:
K_C=K_L+<X,AX>, r_C=r_L-QAX.
For arbitrary response x in Q, set E_C=2 Re<r_C,x>-<x,Cx>, e_C=r_C-Cx. Then the original Schur lower comparison is exactly
    K_C-<r_C,C^-1 r_C>=K_C-E_C-<e_C,C^-1 e_C>.
This uses invertibility of the bounded coercive step comparator C. In particular one must not receive the CR-ALP-069 scalar PASS by merely replacing L^-1 with C^-1 in the old shifted expression.

Let rho be the archived physical proxy to r_L, eps its 058 global L2 error, eta=1e-13 the 067 rational-shift transport for x. Let E_old and y be the saved 067 integrals of 2<rho,x>-<x,(D_old+P_rat)x> and rho-(D_old+P_rat)x. Write e0=Qy. Then
    K_C-E_C >= K_L-E_old + <X+x,A(X+x)>
               - <x,(D_old-gamma)x> - (2 eps+eta)||x||.
This identity retains the useful positive correction energy and its mixed term by forming the square of X+x, not by discarding mixed terms. The signed D_old-gamma multiplier may have either sign; use an outward upper endpoint of its quadratic value. The combined energy error follows from C_old-C_old,rat=Q(P-P_rat)Q and r_L-rho. This direct margin comparison does not require an inverse transport of the possibly large new residual.

## 2. Full residual without new quadrature

For each high ordinary sine coefficient define
    h_j=(D_new,j-D_old,j)x_j + A_j X_j.
Then e_C=Q(e0-h)+u, ||u||<=eps+eta. For ANY real fixed pole subtraction t, weighted Cauchy-Schwarz in the constraint space yields
    <Q(e0-h), C^-1 Q(e0-h)>
         <= sum_(j>85) |e0_j-h_j-t c_j|²/d_j
where c_j=<e_j,b> is the genuine matching pole coefficient. It is unnecessary to project h explicitly, because Q kills the remaining pole vector and e0 already lies in Q. The same inequality holds for complex data with complex t; the archived witness is real. The true inverse penalty is at most (sqrt(U)+(eps+eta)/sqrt(dmin))² if U is a certified upper bound for this weighted sum. This follows from the triangle inequality in the C^-1 norm.

Reconstruct e0 through the 067 saved projected residual formula:
    B=<y,b>-sum_(j<=85)c_j y_j,
    st=||b||²-sum_(j<=85)c_j²,
    e0_j=y_j-B c_j/st for j>85,
    ||e0||²=||y||²-sum_low |y_j|²-|B|²/st.
Then Parseval gives T0=||e0||²-sum_(85<j<=2048)|e0_j|². Use outward upper T0 and finite coefficient intervals.

## 3. Full-pole correction tail

The original 056/057 correction has representation X=sum_(j<=850) D_j^X e_j+Db b; the complete pole function is not truncated. Since x has only ordinary sine coefficients through N=2048, X_j+x_j=Db c_j for j>N. Require original X coefficient enclosures tied to the same 056 rational witness and K_L; archived midpoint proxy coefficients alone are not automatically such enclosures. Check the extraction/assembly source before using them.

For the tail, cTail=sum_(j>N)c_j²=st-sum_(85<j<=N)c_j² exactly. The new diagonal gain term has lower tail bound
    sum_(j>N) A_j |X_j+x_j|² >= min_tail A_j * Db² * cTail.
Keep it, together with the exact finite lower endpoint of sum A_j|X_j+x_j|². The tail of the weighted residual cost is bounded by
    [sqrt(T0/dtail)+sqrt(cTail * max_(bands intersect j>N) |A_band Db+t|²/d_band)]².
This follows by weighted triangle inequality, so does not assume any sign or phase of the unmeasured e0 tail. It retains the full pole contribution, even beyond the last band. For each band the denominator d_band is constant and positive; use the actual max across all tail bands including the infinite last band. Finite contribution U_finite=sum_(85<j<=N) |e0_j-h_j-t c_j|²/d_j. A valid U is U_finite plus that tail bound. Any finite set of rational t choices may be evaluated, retaining the smallest rigorous U; no optimization proof is necessary. One may also refine tail pole masses by exact band sums or choose a new trial if the conservative triangle bound fails. This first receiver is deliberately cheap, no new quadrature. It is a sufficient lower bound only; a nonpositive result is UNDECIDED, not an original negative witness.

## 4. Executable dependency and custody

Original CR-ALP-056/057/058 archive reassembler at src/audits/cr_alp_058_archive authenticates the full original ZIP. New source-only recovery workflow at f877a624bff5f3996589c566aa60b88dc7a11fbe checks stage SHA256SUMS, replays 056 extended and 058 endpoint, prints the actual archived prepare.py and assemble.py, and extracts the original 425 even ordinary sine correction intervals plus Db from input_2.txt. The latter extraction is a candidate input bank until its relationship to the original exact X and the witness contraction is checked against the recovered assembly source. It is not the separate 982-coordinate response bank. No numerical receiver result yet. The independent 069 native audit run 45625224037 was queued at last check. Do not promote any original even direction before the receiver margin is positive and its source linkage verified. Prior odd witness, 1048 positive dimensions, conditional index cap 66 unchanged; full 85-column Schur, arbitrary support, Lean formalization and RH unresolved. Historical standalone archive reconciliation open.
