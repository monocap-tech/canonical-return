# CR-ALP-042 — Support-aware prime translation norm reduces the positive-tail cutoff to 172

**DRAFT / NONPROMOTING.** New classical analytical inequalities with executed exact integer and Fraction arithmetic. No Lean compilation, full-collar positivity, arbitrary-support transfer or RH claim.

## 1. Source contract and physical support

In the original fixed-collar two-pole-null problem, a=11/10, I=(-a,a), and the source is `Zeta23.EF.literatureRHS` (weil-lab `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, pinned blob `0034976fe2f1e5fa9d411108736160542eb94293`). The genuine Gamma symbol is g(r)=Re psi(1/4+ir/2)-log pi. On the two-pole-null space H0 both pole terms vanish. The prime part P is minus the compression to I of the seven positive paired translations with weights omega_n=Lambda(n)/sqrt(n), n=2,3,4,5,7,8,9. CR-ALP-027 bounded ||P|| by 2 sum omega_n<177/25, ignoring the fact that the translated zero-extended sources must remain in the physical interval. We retain precisely the same source operator; the following enlarged interval is used solely for bounding its norm.

Use y=exp(x+a) and the physical measure dx=dy/y. The original interval becomes J=(1,exp(11/5)). Exact atanh-series bounds give log 9<11/5<log 10, so J is contained in U=(1,10). Define the nonnegative weighted translation operator on L²(U,dy/y), with zero extension, by

    (Af)(y) = sum_n w_n [1_U(n y) f(n y) + 1_U(y/n) f(y/n)],

where w_n>=omega_n are rational upper enclosures. Every term in a pair is the adjoint of the other, because dy/y is invariant under multiplication by n. Hence A is bounded and self-adjoint. It preserves nonnegative functions; positive semidefiniteness is neither assumed nor needed. For original sources, pointwise domination by A acting on the absolute zero-extended source implies ||P||<=||A||. This does not replace the source form with an enlarged-collar form or introduce any new prime interactions.

## 2. Exact finite step-cell recursion for the tenth power

The rational upper weights, with denominator 10^6, have numerators

    2:490130, 3:634285, 4:346574, 5:719763,
    7:735485, 8:245065, 9:366205.

They are generated using 40-term rational atanh/log enclosures and integer square-root lower bounds. No floating-point values participate in acceptance.

Let B_0={1,10} be the breakpoints of the constant function v_0=1 on U. Define recursively

    B_k = {1,10} union {n t, t/n: t in B_(k-1), n in powers} intersect (1,10),
    v_k = A v_(k-1).

The endpoint intersection is applied only to the new interior points. Every v_k is constant almost everywhere on each open cell between consecutive breakpoints. This follows by induction: membership in U and the old cell index of n*y or y/n cannot change inside a new cell. Evaluation at the exact rational midpoint therefore determines its entire cell value, not merely a sample estimate. The exceptional endpoints form a finite null set. Computing with integer numerators over (10^6)^k gives 2,189 cells at k=10 and

    ess sup A^10 1 =
    15306858409238293405134184753033878801746910636705348970868988259 / 10^60
    < (21/8)^10.

For completeness, each term of A^10 corresponds to a path of ten paired shifts, with a nonnegative product weight and an indicator that every intermediate coordinate belongs to U. For any f, pointwise weighted Cauchy-Schwarz gives |A^10 f|² <= (A^10 1) A^10(|f|²). Integrating and using self-adjointness (or reversing every path) gives ||A^10 f||² <= C²||f||² for C=ess sup A^10 1. The argument extends from bounded functions to all L² by density. Thus ||A^10||<=C and self-adjointness yields ||A||^10=||A^10||. Consequently

    ||P|| < 21/8 = 2.625.

The code enumerates the exact rational cell partition and checks all row masses. Its comparison C<(21/8)^10 is strict in exact integers. The norm bound is not an estimate of the actual least eigenvalue of P.

## 3. New whole-infinite-tail theorem

Retain the CR-ALP-025/027 physical projected Dirichlet sine space E_N, H0, and Q_N exactly. The classical CR-ALP-027 fifth-order leakage proof gives, for q in Q_N H0 and nu_(N+1)>=2R,

    integral_{-R}^R |qhat|² dr/(2 pi) <= d ||q||²,
    d = 4096 a^5/(45 pi^7 N^5) * (R^5/5 + R^3/6 + R/16).

The CR-ALP-026 digamma bounds give g(r)>-6 globally and g(r)>log R-log(2 pi)-2/R on |r|>=R>=2. Enclose log(2 pi)<37/20, pi>3; round log R down to 10^-6 and d upward to 10^-9, exactly. Write L for the resulting high-frequency Gamma floor and d_up for the leakage upper bound. Then

    W(q,q) >= [L(1-d_up)-6d_up-21/8] ||q||².

All rational comparisons below pass. In particular for N=172, R=117, the certified coefficient is 4927859705274599/19500000000000000 > 1/4. The frequency separation follows from 3(N+1)>=4aR. This covers the whole infinite-dimensional complement in the original logarithmic form domain, not a numerical sample of q.

| N | R | Certified strict lower tail coefficient |
|---:|---:|---:|
| 172 | 117 | 1/4 |
| 200 | 135 | 1/3 |
| 242 | 165 | 3/5 |
| 300 | 200 | 3/4 |
| 15000 | 10000 | 47/10 |

N=172 is an explicit sufficient cutoff, not a minimality claim. The new smaller target is independent of the historical F41 242-mode family: no assertion of equality or transfer is made merely from the dimension.

## 4. Effect on the energy-conditioned CR-ALP-041 comparison

The CR-ALP-041 exact verifier was re-executed with `--verify` on its pinned original signed fourteen-mode interval bank and saved rational basis. Both parity bounds rho<13/125 and rho<79/1000 pass unchanged. At N=15000 the improved tail coefficient 47/10 permits using Young's inequality

    2 |<Tp,q>| <= 5 ||Q_N T p||² + (1/5)||q||².

Writing E for the positive conditioned energy budget and sigma=10^-14 gives

    W(p+q,p+q) >= sigma||p||² + (1-5rho) E + (9/2)||q||².

Both rho bounds are strictly less than 1/5. Therefore the previous mixed E14+Q15000 theorem strengthens to a 9/2 tail margin (strict for every nonzero permitted p+q). No new finite matrix was computed here.

## 5. New acceptance target and limitations

The sufficient full fixed-collar pole-null test now becomes

    H_172 - 4 B_172 >= 0,
    (B_N)_{ij}=<Q_N T v_i,Q_N T v_j>,

in the actual physical basis (or an exactly congruent one). This is UNEVALUATED. The low block beyond the previously certified E14 still requires 158 more finite directions for this target, and all interactions with Q172 must be certified at the new cutoff. The known E14+Q15000 mixed positivity does NOT imply E14+Q172 mixed positivity; positive subspaces cannot be added without checking coupling. Moving the cutoff down may increase residuals and the available trace bound may fail. Alternative N=242 or 300 gives a larger tail margin if that helps a rigorous comparison. No receiving Weil-Lab, F41/F42, formalization, Hodge, provenance or RH acceptance state is promoted.

Next CR-ALP-043: evaluate the original signed finite and residual blocks at the smaller cutoffs, starting with E14 against Q172/Q242/Q300 and the actual high-sine variance. Reuse exact energy-conditioned coordinates, but do not assume CR-ALP-041's old-cutoff residual bounds apply unchanged. The previous 15000 cutoff remains a valid fallback. Full fixed-collar positivity, arbitrary support and RH are open.
