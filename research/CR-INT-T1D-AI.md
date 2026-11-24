# CR-INT-T1D-AI — Pole-null rank-two displacement, direct interval Schur elimination and signed-band mixed-tail control

**Disposition:** DERIVED / ANALYTIC with executed exact-rational interval audit PASS, conditional on AE/AF/AG/AH's original signed-Weil source identities and Q's joint arithmetic norm. No Lean, remote CI, independent mathematical review, global fixed-window positivity, completed Suzuki–Connes transfer, Phase IV or RH.

## 1. Frozen physical source

Work on I=(-2,2), with physical L²(dx), q₀=Γ−A (the original Weil form on the two-Mellin-pole-null source), the full digamma multiplier and all 24 genuine prime powers. Use unnormalized Neumann cosines p_n(x)=cos(nπ(x+2)/4), masses 4 for n=0 and 2 otherwise, parity σ=(-1)^p. In each parity use raw modes n=p,p+2,...,p+82, i.e. 42 raw modes; remove the one pole moment exactly. Write t_i=ω_i², v_i=t_i+1/4, r_i=v₀/v_i and f_i=p_(p+2i)−r_i p_p for i=1,...,41. This gives V82, with 41 dimensions per parity in the physical logarithmic form domain. The two pole terms vanish on the constrained source; no identity off that source for the complete Weil form is claimed.

## 2. Exact rank-two displacement survives pole elimination

Let M be the raw same-parity matrix of the pole-free form with the j=0 boundary term removed as in AG (that term vanishes after the constraint). For i≠j, M_ij=(F_i−F_j)/(t_i−t_j), F_i=F_σ(ω_i), including the entire boundary series j≥1 and all negative arithmetic shifts. Write k_i=M_i0 and d₀=M_00. The genuine constrained Gram H is

    H_ij = M_ij − r_i k_j − r_j k_i + r_i r_j d₀.

For Δ=diag(t_i) on the remaining indices, the exact generators are

    g_i = F_i + v₀ k_i − r_i(F₀+v₀d₀),   h_i=1−r_i,
    ΔH−HΔ = g hᵀ−h gᵀ.

Proof: multiply H_ij by v_i v_j and use v_i−v_j=t_i−t_j. Put u_i=v_i F_i, b_i=v₀v_i k_i, c=v₀²d₀. The scaled H is v_i v_j M_ij−b_i−b_j+c. Its displacement is (u+b)∧v−(v(b−c))∧1, with ∧ the antisymmetric outer product. The relation (v_i−v₀)b_i=v₀u_i−v₀F₀v_i gives v_i(b_i−c)=v₀(u_i+b_i)−(v₀F₀+c)v_i. Thus the displacement equals (u+b−(v₀F₀+c)1)∧(v−v₀1). Divide both generators componentwise by v_i to obtain the displayed identity. It uses exact moment correction, not an approximation to a Cauchy matrix.

Eliminate pivot p=H_kk>0; c_i=H_ik/p. The Schur complement S_ij=H_ij−p c_i c_j has generators g'_i=g_i−c_i g_k, h'_i=h_i−c_i h_k. Indeed (t_i−t_j)p c_i c_j equals (g_i h_k−h_i g_k)c_j−c_i(g_j h_k−h_j g_k), using the two original displacement identities. Update S_ii=H_ii−p c_i². Recover every future off-diagonal from the two generators and distinct t_i. This is exact structured LDLᵀ, O(m²) scalar operations for the factorization, O(m) active generator/diagonal storage (O(m²) if retaining L). No claim that interval bit complexity is quadratic. For fixed moment rows the triangular solves cost O(m²) per row; the conservative physical trace check below solves for all raw coordinate rows and is O(m³). Do not confuse this with a completed matrix-free certificate for all intermediate modes.

## 3. Direct integer-interval factorization passes at 41 per parity

The delivered standard-library implementation stores endpoints as integers on the 10^−220 grid. All arithmetic rounds outward with integer floors/ceilings; all strict comparisons are exact. Trigonometric, Machin, atanh-log and exponential series use explicit Taylor or geometric remainders. The normalized Euler–Maclaurin series of AH is recomputed with K=320, P=80 (through B_160), R=112; t=K+1/4. Gamma tail after R terms is enclosed by rho^(R+1) Z_(2R+3); boundary tail by Z_2/(2t) rho^(R+1)[(R+2)−(R+1)rho]/(1−rho)² plus the positive exponential majorant. Here rho=(ω/2)²/t² and Z_s=t^(s−1)ζ(s,t); the exact EM remainder is bounded by |B_(2P)|(s)_(2P−1)/[(2P)! t^(2P)]. Use each interval g,h,diagonal directly in the structured pivot recurrence, no floating proposal or dense Gram needed. Every pivot lower endpoint is positive in both 41-dimensional spaces. Minimum pivot lower diagnostics: 8.7251783269e−30 (even), 3.3283878135e−29 (odd). These pivots are not physical eigenvalues.

For any functional vector v in the exact constrained basis, forward solve L z=v with directed intervals; v*H^−1 v≤Σ |z_i|_upper²/p_i,lower. Use each raw-coordinate functional (−r for the base coefficient; the other coordinate rows are unit rows) to bound tr(H^−1 G) by the sum of its dual squares weighted by physical cosine masses. This certifies q₀[f]>10^−86||f||² on V82 (reciprocal trace lower diagnostics 2.0416248781e−86 even, 2.2508608645e−84 odd). The sum of raw-coordinate dual norms also bounds the coefficient ℓ¹ norm relative to sqrt(q₀[f]). A separate initial midpoint structured proposal plus rational congruence passed; the final audit instead uses direct interval generators and no mpmath.

## 4. Mixed tail: exact high pole quotient and signed near-cutoff sums

Keep AH's high cutoff N₀=27400, 53 monotone rational bands d(a)≤D_h(n)=mΓ(ω_n)−203/25 with first floor 3268/10⁶>3/1000. All floors and integral/Abel weights are recomputed. Let Y consist of high Neumann form-domain vectors satisfying both Mellin moments. For each f∈V82, AG's complete divided difference gives raw mixed coefficients

    q₀(f,p_n)=Σ_{k≥0} [F_σ(ω_n) S_k(f)−T_k(f)]/ω_n^(2k+2),
    S_k=Σ_m c_m ω_m^(2k), T_k=Σ_m c_m F_σ(ω_m)ω_m^(2k).

For high vectors in Y, a multiple of the high pole moment row (1/4+ω_n²)^−1 pairs to zero within each parity (the constant normalization is common for all high raw cosines). Hence add T₀(f)/(1/4+ω_n²) to the mixed coefficient; the remaining T₀ term is bounded by |T₀(f)|/(4ω_n⁴), not |T₀(f)|/ω_n². This is an exact quotient of the mixed representer on the pole-null high sector, not an estimate on the full unprojected representer. Normalize high p_n by sqrt(2).

For b_σ(ω)=2Σ_{j≥1}(1−σe^(-4a_j))ω²/(a_j²+ω²), compare the decreasing function ω²/(4x²+ω²) with its integrals from 1/4 and 5/4. The integral on [0,∞) is πω/4. Consequently |b_σ(ω)−πω/2|≤5/2+2e^(-10)/(1−e^(-8))<3. Retain the actual signed phase function a_n=2Σ κ_h sin(nπ log h/4). Then

    F_σ(ω_n)/ω_n² = (π/2+a_n)/ω_n + ε_n/ω_n², |ε_n|<3.

Unlike AH's independent Gamma and prime norm bounds, evaluate the nonnegative square (π/2+a_n)² jointly. On the 12 even-start bands between 27400 and 28000, sum the actual square in each parity by interval complex multiplication of the 24 prime phases. The recurrence is reseeded by binary exponentiation every 32 steps (and at each band start), so no uncertified real or complex evaluation enters the result. Bound the weight by 1/(a²d(a)) on each band, and use ω_n≥(785/1000)n. For n≥28000 use AH's decreasing integral and Abel bounds, including the cross term between π/2 and a_n. All nonzero phase denominators separated from zero by directed intervals; all 24 prime powers included. The resulting leading weighted ℓ² norm upper diagnostics are 0.05550262608113 even and 0.05479855714562 odd (AH independent bound 0.08196542600470). These are upper enclosures, not actual mixed-operator eigenvalues.

Use the new direct structured LDL dual bounds for S_k and T_k, k=0,...,23. The geometric remainder from k≥24 is bounded by 2 C L^48 sqrt(q₀[f])(30/ω_n^49+30L/ω_n^50), where C is the computed raw coefficient ℓ¹/energy bound and L an integer exceeding the largest low ω. This follows from |F_σ(w)|<30w and L/ω_n<1/2. Apply the band-weighted inverse tails and Minkowski, including the improved T₀ term. The final exact checks show the squared relative mixed bound (for the corrected representative on Y) <9/10 in both parities. Printed upper diagnostics: 0.8881874986205489 (even) and 0.8590726339711204 (odd). In contrast, the earlier absolute phase estimate at the same cutoff failed the mixed test on all 82 directions; a 114-dimensional exploratory low test had positive even low-block but did not pass its mixed majorant. These are certificate limitations, not negative original directions.

## 5. Restricted infinite-dimensional positivity and next gate

Conditional on Q and the AE/AG source formulas, q₀[y]≥<D_h y,y> on Y. The mixed bound implies |q₀(f,y)|≤sqrt((9/10)q₀[f]<D_h y,y>). Since sqrt(9/10)<19/20, obtain

    q₀[f+y] > (1/20)(q₀[f]+<D_h y,y>)
             > (10^−86/20)||f||² + (3/20000)||y||²

for f∈V82, y∈Y, not both zero. This holds on the original form completion; by restricting smooth sources whose first 27400 Neumann coefficient projection lies in V82, obtain a smooth-source subspace of codimension at most 27400−82=27318 (13659 per parity). Thus n_−(q₀;D⁰)≤27318, not an actual count. The untouched intermediate modes and their couplings are still open. No genuine negative original witness found. Do not promote this result into the Weil-Lab accepted frontier.

**Reproduce:** `python3 ai_bands.py` with `ai_structured.py` and `ai_tail.py` beside it (Python standard library). `audit_output.txt` records the exact executed tests; `float` conversions only print diagnostics. The analytical source/form identities, Euler–Maclaurin remainders, structured elimination identities, high-inverse bound and infinite-tail arguments are written mathematics subject to independent review, not mechanically proved by the audit.

**Next AJ:** exploit the generator update to eliminate intermediate blocks while carrying a posteriori interval widths and mixed dual functionals, or tighten the high inverse by retaining the positive Gamma boundary energy. The present 82-dimensional result is evidence for the structured method, not evidence that every omitted block has positive spectrum. Complete fixed-window pole-null positivity, remaining pole directions, support-independent Weil positivity, Suzuki–Connes transport, Phase IV and RH remain open.