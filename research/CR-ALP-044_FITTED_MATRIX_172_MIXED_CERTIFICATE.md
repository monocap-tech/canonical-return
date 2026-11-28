# CR-ALP-044 — Fitted residual matrix closes E14 plus the entire Q172 tail

**DRAFT / NONPROMOTING.** Classical analytical proof conditional on the inherited continuum identities and pinned finite interval inputs. New exact local arithmetic PASS. Not Lean, not full E172 positivity, not full collar positivity, not arbitrary-support positivity, not RH.

## Source and dependencies

Same original signed form as CR-ALP-043: a=11/10, physical L² metric, genuine two-pole-null space H0, seven negative paired translations at 2,3,4,5,7,8,9 and complete original digamma Gamma term. Source `Zeta23.EF.literatureRHS`, pinned Weil-Lab blob `0034976fe2f1e5fa9d411108736160542eb94293`. Reuse E14 and the exact rank-two corrected high-sine projector Q_N. CR-ALP-041's stored basis/energy audit and CR-ALP-042's support-aware prime/positive-tail audit were rerun and passed. CR-ALP-043 coefficient identities, interval primitives and CR-ALP-038/039 error envelopes remain analytical premises; their proofs are not independently formalized here. The CR-ALP-043 executable used as a module has SHA256 `a56ef90efd6fdd8d23b532def2bde92b300406e3598f9b827f5cab1145416de0`.

## 1. Fitted matrix bound instead of trace

In each seven-dimensional parity sector let v_i be the saved CR-ALP-041 rational conditioned basis, d_i its certified energy budgets, and phi_i=v_i/sqrt(d_i). Then W(p,p)>=sigma||p||²+||z||² when p=sum z_i phi_i, sigma=10^-14. No orthonormality in physical L² is assumed. The CR-ALP-043 signed coefficient h_i(j;172) is already centered by a legitimate source-dependent multiple of the matching pole b. Write c_j=<e_j,b>. For any rational vector t, further subtract t_i c_j from each coefficient: this is annihilated by Q172 globally, not merely on a finite band. Let K=30000 and let U have rows (h_1(j),...,h_7(j)) for matching parity j in (172,K]. Let c be the corresponding vector of c_j. Accumulate exact outward intervals for A=U^T U, v=U^T c and s=c^T c. Here s>0. Choose each rational t_i by rounding the midpoint of its interval v_i/s to 12 decimal places. No exact optimality assumption is needed. The fitted Gram is exactly

    F_K=A-tv^T-vt^T+s tt^T.

The bound on the actual projected full residual is obtained from the entire centered sine square sum, not from subtracting a finite-band fitted norm from a full norm. Opposite parities are orthogonal. This procedure applies to complex coefficient vectors by Hermitian extension of real symmetric inequalities.

Reuse CR-ALP-043's 10^36 outward-integer digamma, exponential and trig primitives, with periodic interval rotation restarts, now through 30000. All source signs combined before products. There are 14914 rows in each parity. Certify F_K < b I using a rational LDL on the interval midpoint of bI-F_K minus its maximum absolute row radius times I. The row radius bounds the norm of every admissible symmetric interval perturbation. This passes with b=9/50 (odd sine indices) and b=4/25 (even sine indices). Independent saved-matrix checker uses exact integer Bareiss leading determinants instead of the primary rational LDL. Neither checker assumes residual columns are orthogonal.

## 2. Infinite remainder: preserve fitted subtraction and endpoint term

For j>K, express each normalized phi_i as an ordinary sine combination sum_l a_il e_l minus g_i b, where g_i is the normalized anchor coordinate times its actual gamma_k. Denote w_n=Lambda(n)/sqrt(n). From CR-ALP-043 and (39.2)-(39.3), after all legal centerings the pole coefficient equals

    C_b(j) = c_j [ (1/2) log(j/N) - sum_n w_n cos(nu_j log n) ] + epsilon_b(j),
    |epsilon_b(j)| <= C2/j²+C3/j³,
    C2=(121/225)[(12001/4000)B0+(11/7)B1]+2 B1 S_a(2a/3)^2,
    C3=8a^3 B0/945.

Here N=172, S_a=177/50, B0,B1=(6/5,3/10) or (3/5,3/5). The extra C2 term bounds the genuine prime sine remainder in the exact translated-pole formula of CR-ALP-043; sqrt(a)>1 and nu_j>3j/(2a). The CR-ALP-039 centered-log error estimate is uniform for j>15000; changing N only changes the explicit log term, so it remains valid here. The ordinary sine combination has the CR-ALP-038/043 signed-moment error envelope A_i/j²+B_i/j³. In the notation D=sum a_il nu_l and C_s=sum |a_il| nu_l^s (s=0,2,3),

    A_i=(2 S_a/a+(201/200)(22/7)/(2a))(2a/3)^2 |D|,
    B_i=(S_a C_0/(2a)+8 S_a C_2/a+(201/200) C_3/a)(2a/3)^3.

The moment symbols C_s here refer to ordinary sine moments, not the pole error constants. The same proof covers a non-pole-null ordinary sine component because the centered subtraction removes the beta=1/2 term. The exact signed leading pole contribution after the fit is c_j[-g_i/2 log(j/N)-t_i+g_i sum_n w_n cos(nu_j log n)]. Let C=2|b(-a)|/(sqrt(a) (pi/(2a))) so |c_j|<=C/j. Write A'_i=-g_i/2 and B'_i=A'_i log(K/N)-t_i; this preserves the cancellation at K before absolute values.

For m=0,1,2, the parity grid moment bound is

    sum_{j>K, matching parity} log^m(j/K)/j² <= J_m=m!/(2K)+2/K².

Proof: the nonnegative function on [K,infinity) is monotone or unimodal, with integral m!/K and supremum <=1/K². On each monotone side of the maximum, compare the mesh-two terms with disjoint adjacent intervals of length two on that side; at most two terms around the maximum/boundary are left over. This gives the stated bound. The two parities need not have the same first node. Thus the squared Frobenius norm of the logarithm-plus-fit part of the remainder is bounded by

    L²=C² sum_i [(|B'_i|² J_0)+2|A'_i B'_i| J_1+|A'_i|² J_2].

The oscillatory prime part has squared Frobenius norm <= P²=C² S_a² J_0 sum_i |g_i|². This uses the full seven-weight sum and does not discard oscillatory jumps. The faster remainder has squared Frobenius norm <= E²=sum_i [(A_i+|g_i| C2)/sqrt(3K³)+(B_i+|g_i| C3)/sqrt(5K^5)]². The code rounds square-root denominators down, numerator intervals outward, and all norms upward. By the triangle inequality for the Frobenius norm (which upper bounds operator norm), the full infinite remainder Gram is <= (L+P+E)² I. This is conservative; no oscillatory orthogonality assumption is needed. The disjoint finite and infinite bands add in squared norm, not in norm.

## 3. Executed results and consequence

| sine parity | finite Gram upper | infinite remainder upper | total rho upper | available mu |
|---|---:|---:|---:|---:|
| odd indices | 9/50 | 30537/1000000 | 210537/1000000 | 1/4 |
| even indices | 4/25 | 219/8000 | 1499/8000 | 1/4 |

Both exact totals are < eta=23/100 < mu=1/4. Initial K=15000 run failed the sufficient odd comparison (total upper 54547/200000>1/4); preserved failure log. Increasing K to 30000, rather than changing original source, resolved it. The saved interval Gram verification also passes all seven Bareiss determinants per parity, using independently rounded component norm budgets to recheck the final strict inequality. It does not independently rebuild all coefficient enclosures or the underlying continuum identities.

For p in E14 and q in Q172 H0 intersect D_log, Young's inequality with eta=23/100 gives

    W(p+q,p+q)>=10^-14||p||²+(1-rho/eta)||z||²+(1/4-eta)||q||².

Thus for nonzero p+q,

    W(p+q,p+q)>10^-14||p||²+(1/50)||q||².

The mixed positive subspace has codimension 172-14=158 in H0 (in the form domain). Since E172 minus E14 consists of exactly 158 orthogonal physical shell directions, the nonpositive index (dimension of any linear subspace on which W<=0) is <=158, conditional on the continuum inputs, by the dimension/intersection argument. This improves the Q172-alone bound 172. It does not prove positivity in those missing directions, or even their isolated finite section. No relation to the frozen F41 basis follows from dimensions alone.

## 4. Reproduce and next

Run from the package directory:

    python cr_alp_041_energy_conditioned_schur_audit.py --verify
    python cr_alp_042_support_aware_prime_tail_audit.py
    python cr_alp_044_fitted_matrix_audit.py
    python cr_alp_044_check_result.py

All four were executed locally with passing final results. Only standard-library Python is needed. The new sweep reuses the predecessor's coefficient interval generator and its inherited analytical bounds, not a new Lean or independently reconstructed source formalization. New exact result stores all 98 fitted interval Gram entries and rational fits. No acceptance-ledger change. Next CR-ALP-045: use the 172 cutoff as a working target, construct a larger actual signed low block with near-null-aware arithmetic; use frequency-dependent inverse coercivity only if rigorously transferred to the identical source space. The remaining 158 directions and all their mixed interactions must still be certified.
