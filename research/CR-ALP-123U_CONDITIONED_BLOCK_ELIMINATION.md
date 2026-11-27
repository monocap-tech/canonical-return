# CR-ALP-123U — Conditioned finite block elimination after independent correction

123S run 45927232333 is executing; 123T is staged for its signed output. No new physical evaluation in this note.

## Projected bank conditioning

For each parity let h=(I-S85)b, st=||h||², c_j=<b,e_j>. In the E850 auxiliary bank the vectors are Q85 e_j=e_j-c_j h/st, so their Gram is G=I-cc*/st. Its least eigenvalue is exactly the omitted pole mass fraction:

 lambda_min(G) = (sum_{j>850, parity} c_j²)/(sum_{j>85, parity} c_j²).

Using c_j=C_b v_j/(v_j²+1/4), 15j/11 <= v_j <= 10j/7, and j>85, (1+1/(4v_j²))^-2 >= 99/100. For a same-parity first index m, decreasing integral comparisons give 1/(2m) <= sum_{r>=0} (m+2r)^-2 <= 1/(2(m-2)). Use m<=852 in numerator and m-2>=84 in denominator. Thus the ratio is at least (99/100)(7/10)^2(15/11)^2(84/852)>2/25; last inequality checked by exact integer cross multiplication. C_b cancels; it is nonzero in each parity.

The physical operator C has ell=51/200 and u=7.072938456. Its compressed coefficient matrix C_B therefore satisfies (51/2500) I <= C_B <= u I, so condition number <347 (exact rational comparison). Distinguish this bound from the physical operator condition <28. This justifies well-conditioned numerical proposal generation, not use of its unverified numerical solution as evidence.

## Why a fixed-source auxiliary correction can leave a large finite covector

In original V-conditioned finite coordinates and the projected high bank, write the capped joint quadratic matrix as M=[[H,R*],[R,C_B]], residual vector at the existing normalized trial as (s,e). Here R maps low coordinates to high covectors. Its low source coordinate is allowed to vary within the pole-null span without changing normalized pole moment.

For correction (z,v), energy change exactly
  2 Re(s*z+e*v)+z*H z+2 Re(v*R z)+v*C_B v.
Eliminate v=-C_B^-1(e+R z). Then
  E_new=E-e*C_B^-1 e + 2 Re[(s-R*C_B^-1 e)* z]+z*L z,
  L=H-R*C_B^-1 R.
By the conditional full pole-null comparator positivity, L>0 in finite coordinates. Solve for z using this residual covector and L, then recover v. This is an exact block identity (no assertion that computed numerical Schur complement is accurate).

The independent 123S correction corresponds to z=0, v=-w. Even if e'=0 on the bank, its low residual becomes s-R* w. Thus a larger tail bank alone need not eliminate finite coupling. The next joint refinement can use the bounded C_B solve with 43/42 low RHS columns plus the residual RHS, followed by a much smaller conditioned Schur solve rather than the original nearly-null unconditioned joint matrix. All coordinates must be rounded and original outward intervals used to enclose energy and both residuals. A negative lower comparator energy is not original Weil negativity.

Implementation requires the actual original V-conditioned low and mixed matrices (not replace H by the certified corrected-energy lower K). Assemble mixed columns by original Hss/Hb complete bridge, then apply V with enough precision to retain cancellations. NumPy binary64 may be adequate for the high block proposals but not the V transformation of unconditioned columns; form the conditioned columns in interval or high-precision arithmetic before converting to numerical proposals. Full candidate evaluation uses original coefficients and exact rational trial, not numerical block products.

If 123S is inconclusive, use its actual signed response and energy readback to determine whether the remaining bottleneck is the finite covector or the complete high tail. No complete pole sign or global proof claimed.
