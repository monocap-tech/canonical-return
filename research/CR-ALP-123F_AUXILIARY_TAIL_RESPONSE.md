# CR-ALP-123F — Certified joint response by auxiliary tail trial

Research-only, conditional on the actual common capped comparator and the 122P/R attachment. 123E complete-bridge run 45920542745 is executing. No new runner activated in this note.

## 1. The actual receiving inequality, not a global residual norm

In a fixed parity let V, X, R, C, K, B_C and T be precisely the 122B/L/P coordinate maps and comparator. Write S_C=K_C-B_C, r_C=R-CX, M=C^-1 R=X+C^-1 r_C. The terminal weighted compression certificate implies

    S_C >= K-T^-1 >= m I, B_C <= T^-1 <= K <= kappa I.

The joint pole-null lower comparator has completed square

    L(z,q)= <z,S_C z> + ||C^(1/2)(q+Mz)||^2.

For a normalized actual pole trial f with finite sine plus full pole coefficient, let a=V*W_phys f, r=Q85 W_phys f. The pole contribution annihilates these test vectors. The response energy for the actual W is bounded by the inverse quadratic form of L, which is exactly

    <a-M*r,S_C^-1(a-M*r)> + <r,C^-1 r>.

This identity follows by an invertible triangular substitution on the physical product Hilbert space. The original supremum is over the log form domain, so enlarging it to that product space is legitimate. No assumption of smoothness of C^-1 r.

## 2. Avoid the worst-case M-norm amplification

Choose ANY auxiliary tail trial y in Q85, and let e=r-Cy. Define

    E_y=2 Re<r,y>-<y,Cy>
    b_y=a-X*r-r_C* y.

Then exact identities give <r,C^-1 r>=E_y+<e,C^-1 e>, and a-M*r=b_y-r_C*C^-1 e. If U is any certified upper bound on <e,C^-1 e>, then

    ||r_C*C^-1 e||^2 <= ||B_C|| <e,C^-1 e> <= kappa U,

by Cauchy-Schwarz in the C^-1 metric. Consequently an explicit sufficient condition is

    sigma >= W(f,f) - E_y - U - (||b_y||+sqrt(kappa U))^2/m.

This differs from 123C: the difficult cancellation is evaluated in b_y (finite signed bilinear products), while only the auxiliary solve ERROR e incurs amplification. It does not charge ||M||^2 against the entire original residual r. A rational Young split with eta>0 yields

    sigma >= W(f,f)-E_y-(1+eta)||b_y||^2/m
                      -[1+(1+1/eta)kappa/m] U.

The eta=1 coefficients using U<=||e||^2/c, c=51/200, are bounded upward by:
parity 1: 2/m <=3973163, 1/c+2kappa/(mc)<=30952827;
parity 2: 2/m <=424018, 1/c+2kappa/(mc)<=3525144.
The rational components 1/m, 1/c, kappa/(mc), 2kappa/(mc) were independently evaluated from the durable 122P result via JS BigInt cross multiplication; displayed combined bounds add the two separate ceilings. They are sufficient coefficients for the auxiliary residual, not revised global physical coercivity.

Even sharper: use S0=K-T^-1 as a rational SPD lower matrix and G0=T^-1 as upper response Gram. Then the transformed finite response is at most (||S0^-1/2 b_y||+sqrt(||S0^-1/2 G0 S0^-1/2|| U))^2. The finite matrix norm factor can be certified by a separate exact generalized Loewner test G0<=lambda S0. No need to construct a floating inverse of the actual physical C.

## 3. Frequency-dependent upper bound with the full pole constraint

The capped tail comparator obeys C>=Q diag(d_j) Q, d_j=Dcap_j-21/8>0, j>85 matching parity. Let e_j^raw be the sine coefficients of ANY representative whose Q85 projection is e (discard low sine coefficients). Let c_j=<e_j,b> be the exact pole coefficients. Weighted completion of the square on the pole-null constraint gives

    <e,C^-1 e> <= U_diag
      = sum |e_j^raw|^2/d_j
         - |sum conj(c_j)e_j^raw/d_j|^2 / sum |c_j|^2/d_j.

All sums converge since d_j>=c and e^raw,c in l2. Proof: maximize 2 Re<e^raw,v>-sum d_j|v_j|^2 over sum conj(c_j)v_j=0; the diagonal inverse solution is projected by the weighted pole functional. Formula invariant under addition of any scalar multiple of c_j to e^raw. It uses the complete high-pole mass rather than truncating its denominator. For a computable upper bound, enclose the full weighted squared norm from above, the weighted pole inner-product modulus from below (zero is always valid), and weighted pole mass from above. With the constant cap beyond 2048 the remaining norm and pole-pair contributions can be attached through full physical Parseval identities, as in 091/111. Do not subtract a finite lower bound on the denominator.

The auxiliary trial y can be selected by the same finite Q85-projected sine bank as 123E. To implement the bound, retain exact rational coefficient witnesses, original normalized trial f, y, b_y, E_y, and a complete upper bound on the auxiliary residual. Existing grouped physical source generators would need a new scale-sensitive error budget: their absolute 5e-9 residual approximation allowance is far larger than these pole energies' square roots. This is a substantive blocker, not cured by replay or changing units. Direct signed coefficient summation with new tail estimates or a new high-precision physical evaluation is required. A normalized rescaling of a residual is legitimate only AFTER the correlated cancellation has been enclosed.

## 4. Status

This is an exact classical conditional receiving reduction. It does not establish the sign of either complete pole Schur scalar or the adequacy of the capped comparator at infinite dimension. The 123E finite census is separate. No native Lean compilation, numerical auxiliary solve or new physical enclosure claimed.
