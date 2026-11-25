# CR-INT-T1D-BL — Constraint-isomorphic Neumann lifts and first finite complement correction

**Status:** executed local integer-interval audits PASS. Conditional negative-index bound remains 240 (121 even +119 odd), with BK's complementary coercivity margins 17/10000 and 53/10000. Complete original signed finite lift matrices are positive, as are finite extensions with four complement directions per parity. No infinite Schur lower-bound or full positivity certificate. Conditional on inherited original fixed-window (-2,2), smooth double Mellin-pole-null source identities, AF Neumann representation, BK form/core transfer and prior audit assumptions. No Lean, remote CI, independent mathematical review or new frequency scan.

## 1. Exact lift coordinates

Let F=Pi T N be BK's physical exceptional synthesis, G=F*F. In each parity p, take unnormalized Neumann modes p_n(x)=cos(n pi(x+2)/4), w_n=n pi/4, n0=p. Set v_n=p_n-r_n p_n0, r_n=(1/4+w_n0^2)/(1/4+w_n^2); exact original pole-null conditions follow from AF. Candidate indices n=p+2,p+4,...<8192. Pivoted QR of numerical F*V proposes 121/119 modes. For the positive dictionary index j, its overlap with p_n is (-1)^((n+p)/2)*(4/pi)*[sin(pi(j-2n)/4)/(j-2n)+(-1)^p sin(pi(j+2n)/4)/(j+2n)], using limit pi/4 at zero; multiply by sqrt(b_j). Constant j=0 even uses sqrt(2 b0) when n=0, zero otherwise. N is exact rational grid. Subtract the r_n reference overlap before contraction by N*. Directed 10^9 grid and checked int64 products enclose H=F*V. For rational inverse proposal Z, ||I-ZH||<1 proves H invertible; W=VH^-1 G then has F*W=G exactly. No claim that W is numerically evaluated to the near-null energy precision. All 22 BK digests match; repeat with strengthened overflow guards reproduces both lift NPZ hashes.

| Parity | Dimension | Maximum n | H enclosure error | inverse residual |
|---|---:|---:|---:|---:|
| Even |121|3590|1.091263e-5|6.070101e-5|
| Odd |119|3349|1.471061e-5|8.820993e-5|

Both exact residual bounds <1/10000. Selected indices and actual rational-grid H, Z saved.

## 2. Complete signed energy evaluator

AF's original unnormalized Neumann Gamma form is diagonal mass*mGamma(w) plus sum of 2 a_k^2(1-sigma exp(-4a_k))/[(a_k^2+w_m^2)(a_k^2+w_n^2)], a_k=2k+1/2, sigma=(-1)^p. For s(w)=sum_k 1/(a_k^2+w^2), use s(w)=Im psi(1/4+iw/2)/(2w) for w>0 and s(0)=psi1(1/4)/4. Full nonexponential boundary term is 2*(w_m^2 s(w_m)-w_n^2 s(w_n))/(w_m^2-w_n^2) off diagonal, and s(w)+Re psi1(1/4+iw/2)/4 on the diagonal (at w=0 this equals 2s(0)). These identities follow from the convergent digamma series and differentiation. Include all 24 signed prime powers using AF's exact product-to-sum matrix entries. Apply exact pole-null congruence. The k=0 exponential term vanishes identically; for k=1..39 retain the signed exponential rank-one correction and bound each omitted entry by 8 exp(-322)/(1-exp(-8)), since projected reciprocal factors are bounded by 2/a_k^2 and a_k>=1. Thus no truncation of the nonexponential boundary series.

Use integer-directed intervals on 10^-110 grid. Evaluate complex digamma and real trigamma by shifting z to Re z=96+1/4, Euler–Maclaurin through B64, and exact recurrence. Remainder bounds for psi and psi1 respectively |B64|/[64 (Re z)^64] and |B64|/(Re z)^65, from periodic Bernoulli bound and |z+t|>=Re z+t. Machin arctan, range reduction, rational Taylor log/atan/trig/exp with explicit tails; no numerical special-function values accepted. Independent AF 2x2 interval comparison passes four entries in each parity. Max entry widths 3.453239e-73 even and 3.481110e-73 odd. Actual complete interval matrices stored as decimal integer endpoints.

## 3. Cancellation-resistant finite sign certificate

Initial 10^-65 run: direct interval LDL stopped after 12/13 positive pivots; exact congruence using midpoint preconditioner also inconclusive due uncertainty. Final 10^-110 direct LDL still stopped after 20 in both parities. Neither is negative evidence. Propose fixed-point midpoint LDL, triangular inverse and diagonal scaling; verify Z is lower triangular with positive diagonal. For integer Z and interval midpoint M and radius E on common grid D, compute Z M Z* and |Z| E |Z|* exactly with arbitrary-precision integer products; directed Gershgorin gives the enclosed congruence after division by D^3. All 121/119 midpoint pivots positive, but only congruence accepted. Both exact bounds satisfy 999/1000 I < Z(V*LV)Z* <1001/1000 I (scaled Z). Preconditioner uncertainty max row <2.52e-16 even /4.71e-15 odd. These establish positivity of the selected finite original-energy lift spaces, NOT positivity of the exact Schur complement.

## 4. First finite complement correction

Add four previously unused Neumann indices 200,400,800,1600 (even), 201,401,801,1601 (odd), heuristic pilot. Let V_n synthesize these pole-null modes, H_n=F*V_n. Exact U=V_n-V C, C=H^-1 H_n, lies in BK's complement, has four independent columns by independence of the selected Neumann modes and reference. H_n enclosed; existing H and energy intervals unchanged entrywise. Complete extended signed energy matrices of dimensions 125 and 123 independently pass the same exact congruence (999/1000,1001/1000). Total 248 finite positive directions, with four complement directions per parity. If A=V*LV, B=V*L V_n, T=V_n*L V_n, define K=B*-C*A and D_U=T-B*C-C*B+C*A C (stars denote adjoints). Exact finite correction W_V=V-U D_U^-1 K preserves F*V. Its energy S_fin=A-K* D_U^-1 K is positive by the extended finite compression. In V constraint coordinates, the full Schur matrix is <= S_fin; in BK's F coordinates congruence by H^-1 G transfers both sides. The finite matrix does NOT bound the full Schur from below.

Propose rounded C0; certify ||C-C0||_F <= ||Z||_F/(1-rho)*(eta_Hn+eta_H ||C0||_F+||H0 C0-Hn0||_F), where rho is the certified inverse residual. Propagate Frobenius errors into K and D_U via exact grid products; Gershgorin bounds D_U, and trace(K*D_U^-1 K) >= (||K0||_F-eta_K)^2/lambda_upper(D_U). Executed exact checks:

| Parity | C error | K error | D_U error | D_U lower | D_U upper | finite gain trace lower |
|---|---:|---:|---:|---:|---:|---:|
| Even |7.244134e-5|.004421|.003312|8.478341|12.485237|2.711329|
| Odd |.000119036|.007928|.007418|4.925727|9.485279|.989519|

Exact conservative bounds 8.4<D_U_even<12.5 and 4.9<D_U_odd<9.5 in operator order; trace gains >2.7 and >.98 in unnormalized V coefficient units, not physical Rayleigh bounds or counts. No numerical values of exact corrected lift claimed beyond these enclosures; exact response is defined by the true matrices. Finite gain positive, full complement residual unevaluated. No inertia subtraction.

## 5. Reproduction and next BM

With bundled BK input checkpoint and AF evaluator in sibling directories, run in cr_int_t1d_bl, for p=0,1: OPENBLAS_NUM_THREADS=4 python bl_lifts.py p; python bl_energy.py p; python bl_precondition.py p; python bl_lifts.py p --extend; python bl_energy.py p --extend; python bl_precondition.py p --extend; python bl_corrector.py p. Then python bl_crosscheck.py. Initial failed audit logs and initial evaluator retained; final files and logs in SHA256 manifest. Full inherited BK proof still requires its upstream dependencies; BL runs from its actual compact checkpoint. Carrier modes selected heuristically, not optimal. Local arithmetic only, no independent formal verification.

**Next BM:** enclose the original full operator action on the constraint-preserving corrected lifts, including the infinite Neumann residual tail and the physical projection P_C. For true residual r=P_C L W, use BK's bound S>=W*LW-mu^-1 r*r or the retained positive-metric inverse refinement. Work in a relative preconditioned energy metric, since raw interval LDL and coarse absolute errors have already failed. Expand finite corrector space only where needed. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.