# CR-INT-T1D-AK — Joint energy-dual Gram majorant

**Disposition:** analytical refinement with executed exact-integer interval audit, conditional on AE–AJ source identities, tail estimates and Q's joint arithmetic norm. Not independent review, Lean, remote CI, complete fixed-window positivity or RH.

## 1. Preserve AJ's comparator and source

Use I=(-2,2), exact Mellin pole-null constraints, the complete Gamma multiplier and all 24 original signed prime powers. Keep N0=27400, the 53 high frequency bands, and AJ's infinite boundary retention j>=8000. The high comparator energy is e[y]=d[y]+b[y]. AJ supplies a separated expansion of the mixed functional in the direct-sum high energy dual space, with 48 terms (S_k,T_k, 0<=k<24) and a remainder. Its term vectors on the high side have certified norms bounded by weights a_i (AJ's leading residual plus boundary dual cost, other moment tails and the pole-quotient T0 bound). Its real low moment rows are v_i. The old upper bound on the mixed norm is sum_i a_i ||v_i||_{H^{-1}} + R, where R is the unchanged geometric remainder bound. No new high-frequency identity or boundary coordinate is assumed in AK.

## 2. Joint Gram bound

For the exact low energy matrix H=L diag(p) L^T, direct AI interval generator elimination provides lower pivots d_i<=p_i and interval enclosures of L. Define z_i=diag(d)^(-1/2)L^(-1)v_i. This is the dual metric of the positive lower comparator L diag(d) L^T <= H. It is legitimate to bound the mixed operator in this stronger dual norm, even though individual off-diagonal entries of the actual inverse Gram need not be dominated entrywise. Define the signed Gram

    G_ij = a_i a_j sum_k (L^(-1)v_i)_k (L^(-1)v_j)_k / d_k.

The finite mixed operator, after whitening by the lower comparator, has Hilbert–Schmidt norm squared at most sum_ij |G_ij|. Indeed, writing its separated terms as z_i tensor h_i with ||h_i||<=a_i, its Hilbert–Schmidt square is the sum of products of the signed low and high Gram entries. Cauchy–Schwarz on the high factors bounds each high Gram entry by a_i a_j. Equivalently the same bound follows from expansion of the finite-rank operator's trace. This argument works in the complexification, since the low rows are real. The operator norm is no larger than the Hilbert–Schmidt norm. Adding the unchanged remainder gives the exact sufficient bound

    beta <= (sqrt(sum_ij |G_ij|) + R)^2.

By Cauchy–Schwarz for G this is no worse than the corresponding separate-moment triangle bound with exact same comparator, before interval overestimates. Unlike the old triangle estimate, it exploits cancellation in the signed low Gram inner products before absolute values. It does not evaluate signed high cross-correlations, or optimize the boundary ansatz; those remain distinct future improvements.

## 3. Exact arithmetic execution

Use the unchanged AJ high_precision dependency modules (10^-320 integer grid). For each parity, direct structured LDL encloses L and all positive pivots. The forward solves are interval-valued. Scale each coordinate by its exact rational norm bound before forming the Gram; do not prematurely round a tiny weight before multiplying large moments. Products and divisions by the rational lower pivots are outward interval operations. Sum the signed coordinate intervals before taking a maximum absolute endpoint. The upper square root uses AJ's rational integer-square-root routine. The entire sign assertion is exact; printed diagnostics are floats only.

Raw parity size 48, hence 47 pole-null dimensions each: all 94 pivots positive. Even joint upper diagnostic 0.9767318562400629 (old triangle 1.0086785802619227); odd joint 0.9559607154430081 (old triangle 0.9872933739823821). The interval physical trace reciprocal diagnostics are 3.7365474065428045e-95 and 6.200211017789437e-93. Exact assertions in the successful rerun give beta<49/50 and physical mass trace <10^96. The even old triangle upper diagnostic exceeding one is not a negative direction, but the joint certificate succeeds where that majorant does not.

Raw parity size 50, hence 49 pole-null dimensions each, has positive interval pivots but joint upper diagnostics 1.0409909212280908 and 1.0225993319981876, compared with triangle 1.0790542794141618 and 1.0602006742729198. The actual sufficient gate failed in both exploratory runs; no negative original direction follows.

## 4. Source consequence

On V94, q0[f]>10^-96 ||f||². The high sector Y remains exactly pole-null with frequencies n>=27400. |q0(f,y)|² <= (49/50) q0[f] e[y]. Since (99/100)^2>49/50, completing the square gives

    q0[f+y] >= (q0[f]+e[y])/100.

Thus conditional strict positivity holds on V94 plus Y. Conservative codimension as in AJ gives negative index <=27400-94=27306. This does not control intermediate modes or their mixed interactions, and no global sign, all-support theorem or Suzuki–Connes transfer follows.

## 5. Reproduction and next gate

Run `python3 ak_joint.py 48 0` and `python3 ak_joint.py 48 1`. Exploratory `50` runs deliberately fail the final sign assertion. All prior analytic estimates are imported, not mechanically proved by these scripts. The bundle includes code and logs, with dependencies copied unchanged from AJ. AK's modest dimension gain is evidence of a tighter upper bound, not of convergence to the true full spectral sign.

**Next AL:** compute the actual high-side cross Gram (at least leading and first few moment terms) using signed finite-band products and a rigorous Abel/integral remainder. Alternatively tighten the joint prime-operator comparator or perform block-local high inverses. The current absolute-entry majorant still fails at 49 dimensions per parity; do not extrapolate it as a proof. Independent review, complete intermediate positivity, completed Suzuki–Connes transfer, Phase IV and RH remain open.