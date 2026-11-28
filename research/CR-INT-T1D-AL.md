# CR-INT-T1D-AL — Signed high-side Gram and extended finite-band evaluation

**Disposition:** derived analytical refinement, with executed local integer-interval audit. Conditional on AE–AK's original signed source identities and Q's joint prime norm. Not independently reviewed, Lean-certified, remote CI, full fixed-window positivity, Suzuki–Connes transfer, Phase IV or RH.

## 1. Source and metric

Keep I=(-2,2), both exact Mellin pole constraints, all 24 signed prime powers, and the complete Gamma series. In one parity use AJ's boundary-retaining high lower comparator with N0=27400 and j>=8000. For AL take the diagonal lower comparator explicitly to be the stepwise diagonal with the 53 certified rational band floors d(n), rather than the actual higher Gamma diagonal. On the high pole-null form sector q0[y]>=e[y]=sum_n d(n)|y_n|²+sum_j |<u_j,y>|². This preserves AJ's direct-sum dual representation, whose norm bounds the mixed form relative to e[y]. Use AI's exact pole-null basis and structured interval LDL. The low dual comparator is L diag(lower p) L^T <= H; no elementwise inverse domination is assumed.

## 2. Signed Gram identity and finite coefficients

The separated mixed expansion has S_k=sum_m c_m w_m^(2k), T_k=sum_m c_m F_sigma(w_m) w_m^(2k). Its high coordinates, before dividing by sqrt(2d(n)), are F_sigma(w_n)/w_n^(2k+2) and -1/w_n^(2k+2), respectively. For S0 retain AJ's boundary cancellation, giving (g_sigma(w)+a(w))/w where a(w)=2 sum_h kappa_h sin(w log h), g_sigma(w)=2w sum_{1<=j<8000}(1-sigma exp(-4a_j))/(a_j²+w²). For T0 add the exact null high Mellin row to obtain -1/[4w²(w²+1/4)]. Boundary dual coordinates are present only for S0; their squared norm is bounded by AJ's cost C=(10001/10000)/(2(8000-3/4)). This is an upper bound on the mixed dual cost, not the actual high inverse.

For the first six high vectors (S0,T0,S1,T1,S2,T2), put P0=(785/1000)N0 and scale the high vectors by P0^p for p=(1,4,3,4,5,6), reciprocally scaling the corresponding low rows. For u=P0/w the scaled unnormalized coefficients are

    ( (g+a)u, -u^4/[4(1+1/(4w²))], (b/w+a)u³, -u⁴,
      (b/w+a)u⁵, -u⁶ ).

The finite Gram is the sum of their outer products weighted by 1/(2d(n)), in the appropriate parity, for 27400<=n<48000. Retain all signed products before interval accumulation. Enclose g on each band by AJ's monotone arctangent endpoint bounds, and b(w)/w by pi/2 +/- 3/w (the AI decreasing integral estimate, including the exponential tail). Phase coefficients are computed by directed interval complex recurrence, reseeded every 32 modes. The high grid is 10^-36, converted outward from AI's 10^-320 enclosures; all arithmetic is integer/Fraction. No approximate phase or quadrature value is used as a bound. Add [0,C P0²] to the S0 diagonal for the boundary cost.

## 3. Rigorous infinite Gram tail

At CUT=48000, for a positive nonincreasing weight n^-2/d(n), the absolute cosine-phase sum is <=1/[CUT² d(CUT)|sin(pi s/4)|] for nonzero separated phases s. The same estimate applies with the monotone residual g multiplying the weight, times g(CUT). For phase zero use the bandwise decreasing integral bound. All nonzero prime-phase denominators are checked away from zero. Expand the squared residual amplitude g+a as in AJ; this yields an upper norm r for S0 beyond CUT. For the complete amplitude use pi/2+a with the same bound and add the norm of the <=3/w discrepancy; obtain a norm s for F/w² beyond CUT. For S_k use s/(W CUT)^(2k), W=785/1000; for T_k use the corresponding bandwise integral of w^-2p (T0 has its additional 1/4 factor). Apply the exact P0 scaling. Cauchy–Schwarz gives signed tail Gram intervals [-t_i t_j,+t_i t_j] for i!=j, and [0,t_i²] for the diagonal. The finite high Gram plus these intervals and boundary diagonal encloses the true direct-sum high Gram. The remaining 42 high vectors retain AK's old norm bounds; the k>=24 geometric remainder remains unchanged.

For the low Gram use interval forward substitution and division by certified lower pivots, as in AK. If G_L and G_H are the exact Gram matrices in the respective energy metrics, the separated operator has Hilbert–Schmidt square tr(G_L G_H). Use the interval product upper endpoint for each of the first six signed high entries with its signed low entry. For the other entries use |G_Lij| times the product of high norms. This bounds the trace from above, even for negative cross terms. Add the geometric remainder in norm, then square. The resulting beta is a valid relative mixed upper bound, not the largest eigenvalue of the original form or of its exact Schur correction.

## 4. Executed comparison

At raw dimension 50 (49 exact pole-null directions per parity), both direct structured LDL factorizations pass. The leading high dual norm square upper diagnostics are 0.001310404275673795 (even) and 0.0013035506089984765 (odd). The residual infinite S0 norm square upper diagnostic is 0.00012566263046663535 in each parity. The high Gram resolves definite signs, including negative S0,T0 and positive S0,S1 cross entries. The joint relative bound diagnostics are 0.8656938508325671 and 0.8601305119229297; retaining the same improved norm bounds but taking all high cross entries absolute gives 0.8892247086065893 and 0.8832496106499345. AK at the same dimension had 1.0409909212280908 and 1.0225993319981876. Thus the newly enclosed high norms and signed cross terms together resolve the failed AK gate. This does not establish that the old bound exceeding one reflected a negative direction.

The original two executions assert beta<1 exactly. The reciprocal physical trace diagnostics are 6.063172604469629e-98 and 9.118436643549819e-96; these are lower-bound diagnostics, not actual minimum Rayleigh values.

## 5. Source consequence and next gate

For the 98-dimensional exact pole-null V, all interval pivots are positive. The proven mixed inequality |q0(f,y)|² <= beta q0[f] e[y], beta<1, yields strictly positive q0 on V plus the pole-null high sector. On the smooth core use the inverse image of V under the low Neumann projection as in AG. Conditional negative index <=27400-98=27302. This is an upper bound, not a measured index. Intermediate modes and interactions remain uncontrolled. No original negative direction, whole-window or support-independent positivity, completed Connes/Sonin transport, Phase IV or RH.

Reproduce with `python3 al_signed.py 50 0` and `python3 al_signed.py 50 1`, using bundled AI/AJ dependencies. The local logs record the actual calculations. The mathematical source identities and tail derivations are analytic imports, not mechanically proved by the scripts.

**Next AM:** move beyond incremental increases in low trial size: use a block-local high comparator, an improved joint arithmetic norm or a certified structured factorization with a genuinely reduced intermediate block. The signed Gram improvement should be reused, not extrapolated as a convergence theorem.

## 6. Executed strengthened and enlarged certificates

`python3 rerun.py 0` and `python3 rerun.py 1` recomputed the signed high Gram and the 49-dimensional factorization, checking exactly beta<87/100 and physical trace<10^98. Both assertions passed in both parities. Thus on V98, q0[f]>10^-98||f||², and |q0(f,y)|²<(87/100)q0[f]e[y] for nonzero f,y. Since (15/16)^2>87/100, q0[f+y]>=(q0[f]+e[y])/16.

The same runs then enlarged the raw parity dimension to 54 (53 pole-null vectors each). All 106 interval pivots were positive. The mixed upper diagnostics are 0.9530753849177206 (even) and 0.937850710502224 (odd); the same norm bounds without signed cross entries give 0.9860332726238266 and 0.9690906587528613. Both exact beta<1 assertions passed. Reciprocal physical trace diagnostics are 3.284727398538168e-103 and 1.0048254122032891e-100. No fixed rational physical lower margin or stronger beta threshold was asserted for this enlarged block; the positive pivots and exact beta<1 suffice for positivity of V106 plus Y. The new conditional conservative negative index bound is 27400-106=27294. No further enlargement or full intermediate elimination tested. The signed Gram estimate is a sufficient criterion, not a calculation of the exact mixed operator norm.