# CR-ALP-123C — Direct terminal-coordinate inverse response bound

Research-only. The 123B extended coefficient and physical residual census is executing (45919679528). This analytical reduction uses the existing conditional 122P certificate without another native runner. It does not assert either pole Schur sign.

## 1. Residual functional in the actual receiving coordinates

In each reflection sector write any pole-null test v=Vz+q, where V is the actual conditioned finite coordinate map, q in the Q85 tail form domain. The 122L/P/R lower bound is

    W(v,v) >= m ||z||^2 + c ||q+Mz||^2,

where c=51/200, m=delta lambda/(kappa+delta), M=X+C^-1(R-CX), and ||M||^2<=M2 (saved square_shift_norm2_upper). All are the same physical operators and coordinates as 122P. No claim that the shifted square belongs to the log domain is needed.

For any normalized pole trial f in D with finite sine plus pole representative (hence original signed physical action is L2 by 122A), define a=V* W_phys f and r=Q W_phys f, where W_phys f is the L2 representative of the actual signed form action on this f. For pole-null tests the pole contribution vanishes. With the first-antilinear convention W(f,v)=<a,z>+<r,q>. Substitution y=q+Mz gives W(f,v)=<a-M* r,z>+<r,y>. The full response energy is

    R_f = sup_{v in D0} [2 Re W(f,v)-W(v,v)]
        <= ||a-M* r||^2/m + ||r||^2/c.

Indeed the supremum may be enlarged to all independent z,y in the finite Euclidean and tail physical Hilbert spaces; completing both scalar squares gives the bound. The associated exact Schur scalar is sigma=W(f,f)-R_f since the pole moment of f is one. This formula bounds the response of any actual normalized trial, without requiring f to be a minimizer.

For an arbitrary approximate Galerkin trial,

    R_f <= (2/m)||a||^2 + (1/c+2M2/m)||r||^2.

If the finite Galerkin equations hold exactly, a=0, and

    R_f <= (1/c+M2/m)||r||^2.

The latter is not to be used for a rounded approximate trial without certifying zero finite residual. The former applies to the 123B rational trials, once both the conditioned finite residual and the entire tail norm are enclosed. A finite high-band projected residual lower bound is not such a norm upper bound.

## 2. Exact rational constants

Reopened 122P result blob 7fb998716bdb5a8f0c178585f1f8ad00e7d3fe37 and independently evaluated the following rational inequalities by BigInt cross-products. Rounded upward integer bounds:

| parity index | ceil(2/m) | ceil(1/c+2M2/m) | ceil(1/c+M2/m) |
|---|---:|---:|---:|
| 1 | 3973163 | 11878888495 | 5939444250 |
| 2 | 424018 | 3463241970 | 1731620987 |

Thus admissible simple response bounds are

    cosh: R_f <= 5000000 ||a||^2 + 20000000000 ||r||^2
    sinh: R_f <= 500000 ||a||^2 + 4000000000 ||r||^2

and for an exact Galerkin trial the tail-only coefficients 6000000000 and 1800000000 suffice. In contrast the previous physical coercivity route charges 10^38 and 10^34 times the complete physical residual norm squared. The present gain concerns the tail residual of trials with small finite coordinate residual, not an improvement to global L2 coercivity.

Using the 123A finite Schur energy diagnostics at E85, exact-Galerkin tail norm targets would be approximately 2.65e-24 and 2.08e-21 respectively; these are illustrative, not certified thresholds for the actual rounded trials. Even the sharper receiver may remain insufficient if the true tail residual is larger. A finite upper bound on the Schur scalar never supplies a positive infinite lower bound.

## 3. Next executable gate

After 123B, compare its high-band projected residual lower bounds against the improved tail budget, not only the crude physical coercivity threshold. Independently enclose the finite coordinate residual a using the same 049 V and original 046 matrix/cross intervals. For a valid positive Schur witness require a complete tail norm upper bound or a sharper weighted inverse-response bound, not the finite band alone. If this still fails, investigate tail response with the full capped comparator and a corrected joint finite/pole system. No new native runner, no unconditional collar positivity or RH.
