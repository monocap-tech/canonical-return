# CR-ALP-043 — Exact signed coefficient band closes E14 with Q242 and Q300

**DRAFT / NONPROMOTING.** Classical continuum argument with executed outward-integer interval arithmetic. Not a Lean theorem, not a certificate of the full 172-, 242- or 300-mode finite sections, not complete fixed-collar positivity, not arbitrary-support positivity, not RH. The smaller-cutoff result retains the original source, genuine physical Gram and the two actual pole constraints.

## 1. Source custody and target

Use the fixed physical collar I=(-11/10,11/10), H0=ker B*, e_j=a^(-1/2) sin(nu_j(x+a)), nu_j=j*pi/(2a), a=11/10, and the original signed Gamma-plus-prime operator T from `Zeta23.EF.literatureRHS`, Weil-Lab `Def_Zeta23_ExplicitFormula.lean`, pinned source blob `0034976fe2f1e5fa9d411108736160542eb94293`. Its seven genuine negative paired prime terms have n=2,3,4,5,7,8,9 and weights omega_n=Lambda(n)/sqrt(n). No positive surrogate for the prime form is used.

Use the original projected sine spaces E_N and exact orthogonal projectors Q_N onto their pole-null complements. In each parity sector of E14 use the exact rational coordinate vectors V and positive energies d_i from CR-ALP-041. For p=sum z_i v_i, that theorem gives W(p,p)>=sigma||p||²+sum d_i |z_i|² with sigma=10^-14. Its stored basis/energy audit was rerun successfully against the original CR-ALP-040 matrix and metric inputs. Their SHA256 values are respectively `7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c` and `6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9`.

CR-ALP-042's support-aware prime norm and tail comparisons were also rerun successfully. The analytical lower bounds on W(q,q)/||q||² at N=172,242,300 are respectively mu_N=1/4,3/5,3/4. In this stage the objective is to bound the actual full cross residual, not merely to inherit the old Q15000 bound. No statement about the frozen F41 family follows from the dimension 242.

## 2. Exact original signed coefficients after pole subtraction

Set beta_m=2m+1/2, L=2a, sgn=(-1)^j, b the matching genuine cosh(x/2) or sinh(x/2), b0=b(-a), b1=b'(-a), c_j=<e_j,b>=2 b0 nu_j/[sqrt(a)(nu_j²+1/4)]. For matching parity put

    R(v)=sum_(m>=0) 1/(beta_m²+v²)=Im psi(1/4+i v/2)/(2v),
    R_s(v)=sum_(m>=0) (1-sgn exp(-L beta_m))/(beta_m²+v²),
    Theta_i=sum_n omega_n sin(nu_i log n),
    A_i=nu_i R_s(nu_i)+Theta_i,
    S(v)=sum_n omega_n sin(v log n).

For high j>14, v=nu_j, the complete Gamma resolvent sine cross-entry (CR-ALP-029) and the genuine signed prime translation formula (CR-ALP-038) give the ordinary sine coefficient with a legal pole subtraction:

    C_i(v)=(2/a)[v A_i/(v²+1/4)
                    -(v A_i-nu_i(v R_s(v)+S(v)))/(v²-nu_i²)].

Indeed the second term is the full off-diagonal signed operator coefficient. The first is minus a scalar multiple of c_j. The constant Gamma multiplier contributes zero on the high sine modes of e_i. This formula applies to the ordinary sine component of the anchor as well as to the interior w_i. No pole-null constraint on e_i is used in the formula.

For the pole function itself define

    E_m=int_0^L exp(-beta_m t)b(-a+t)dt-b0/beta_m,
    E_inf=sum_(m>=0) E_m,
    E(v)=sum_(m>=0) E_m/(v²+beta_m²),
    Lambda(v)=Re psi(1/4+i v/2)-psi(1/4),
    t_N=(log(nu_N/2)-psi(1/4))/2-E_inf/b0.

The convergent Gamma coefficient after subtracting t_N b, and dropping the constant g(0)b killed by Q_N, is

    C_b^Gamma(v)=Lambda(v)(c_j-b0/(sqrt(a)v))-2v E(v)/sqrt(a)-t_N c_j.

This is CR-ALP-039 (39.1)-(39.2). It removes the constant logarithmic term rather than deleting any resolvent. The following exact partial-fraction evaluation avoids summing the slowly decaying E_m tail. For m>=1 set s=-sgn and

    E_m=(b0/4+beta_m b1)/(beta_m(beta_m²-1/4)) - F_m,
    F_m=s exp(-L beta_m)(beta_m b0-b1)/(beta_m²-1/4).

For m=0, with ep=(b0+2b1)/2, em=(b0-2b1)/2,

    E_0=ep L+em(1-exp(-L))-2b0.

The elementary sums are

    S1=sum_(m>=1)1/(beta_m²-1/4)=1-log 2,
    S2=sum_(m>=1)1/[beta_m(beta_m²-1/4)]
       =2 psi(5/4)-psi(1)-psi(3/2),
    R1(v)=R(v)-1/(v²+1/4),
    S3(v)=Lambda(v)/(2v²)-2/(v²+1/4),
    E_inf=E_0+b1 S1+(b0/4) S2-sum_(m>=1) F_m,
    E(v)=[E_0+b1(S1-R1(v))+(b0/4)(S2-S3(v))]/(v²+1/4)
         -sum_(m>=1) F_m/(v²+beta_m²).

All nonexponential sums above are complete. The F_m correction and the exponential part of R_s are absolutely convergent geometric-tail sums.

For a genuine translation t=log n, direct integration of sin(vy)b(-a+y+t) over the actual overlap interval, with its reflected partner, gives

    U_jb(t)=2/[sqrt(a)(v²+1/4)] [v(b(-a+t)+b0 cos(vt))-b1 sin(vt)].

The primitive is (b' sin(vy)-v b cos(vy))/(v²+1/4). The constant b(-a+t) term is exactly pole-aligned. Thus the true signed prime coefficient after subtracting that scalar multiple of c_j is

    C_b^P(v)=-2/[sqrt(a)(v²+1/4)] [v b0 sum omega_n cos(v log n)-b1 S(v)].

For each rational conditioned source v_i, express its original physical function as p_i=sum_j d_j^sine e_j-y gamma_k b, where k=1 or 2 is the anchor index and d_k^sine=y-sum_(j>=3) x_j lambda_j. The centered high-sine coefficient is exactly

    h_i(j;N)=sum d_j^sine C_j(nu_high)-y gamma_k(C_b^Gamma+C_b^P).

(The j in C_j is a low source index.) The actual Q_N T p_i norm is no greater than the full square sum of these centered coefficients. We choose the centering t_N explicitly; it need not be the optimal pole variance minimizer. This avoids subtracting uncertain large squares in the exact variance formula. All original signs are combined BEFORE taking interval absolute values.

## 3. Exact finite band calculation

For j>N through K=15000 the new stdlib Python audit encloses the above coefficients in integer intervals scaled by 10^36. There are 7414 high indices per parity in the largest band (N=172). All three cutoffs use the same streamed coefficient calculation. Source coefficients lambda_j, gamma_k and the exponential values are enclosures from the pinned CR-ALP-040 input.

The checker uses rational Machin arctangent for pi, integer square roots, logarithm range reduction followed by 48-term atanh series, and Taylor sine/cosine on a rigorously reduced argument. The rotation recurrence is periodically restarted after 16 indices; interval endpoints propagate outwards. It encloses psi using the classical Euler-Maclaurin asymptotic expansion on the right half-plane, with 12 Bernoulli terms and an explicit remainder. The remainder for Re z>=32, or Im z>=120, is <10^-31: use sup|B_24({t})|<=4*24!/6^24 in the standard integral remainder, and bound int_0^infty |t+z|^-25 dt by 1/(24*32^24) or (1+1/24)/120^24, respectively. Low arguments are shifted by 32 with the exact digamma recurrence. These remainder comparisons are exact rational checks. Each exponential sum is enclosed by 32 terms plus a geometric error <=10^-36, using exp(-a)<2/5. The exact source coefficient widths in the streamed band are below 10^-27. No floating-point value controls acceptance.

For each parity, sum the interval upper bounds of |h_i(j;N)|²/d_i, over all transformed sources and all matching j in (N,K]. Denote the result by rho_band. This is an upper bound on a particular centered residual, not a lower bound from finite-band projection. A separate optional SciPy diagnostic computes finite-band pole variance through index 100000 and is numerical only. It is not used by the exact checker.

## 4. Uniform remainder for every j>K

The same chosen pole subtraction for each ordinary sine component yields the CR-ALP-038 bounds A(p)/j²+B(p)/j³, where D=sum d_j^sine nu_j and C_t=sum |d_j^sine| nu_j^t, t=0,2,3. In particular the proof extends to a sine expansion with an anchor coefficient even without pole-nullity of that sine expansion: subtracting its exact pole-aligned scalar annihilates the m=0 term (beta_0²-1/4=0). Use exactly the CR-ALP-038 formulas for A,B with S_a<177/50. The complete tail norm bound is A/sqrt(3K³)+B/sqrt(5K⁵).

For the b Gamma centered coefficient use the CR-ALP-039 bound (39.3), with the original centering N but summing only j>K. Its log contribution is bounded by

    sum_(j>K) log²(j/N)/j² <= [log²(K/N)+2log(K/N)+2]/K.

The integral comparison holds since K>3N>eN and the integrand decreases there. The remainder terms C2(b)/j²+C3(b)/j³ use the unchanged CR-ALP-039 constants, valid at all j>K=15000. For the exact b prime coefficient above use

    |C_b^P(nu_j)| <= 2 B0 S_a (2a/3)/j + 2 B1 S_a (2a/3)^2/j²,

where B0,B1=(6/5,3/10) or (3/5,3/5) for the two pole functions. The actual endpoint traces and oscillatory translated jumps are retained. Triangle inequality bounds each transformed source's norm on the infinite remainder by r_i. Then rho_tail=sum r_i²/d_i. Both rho_band and rho_tail are rounded upward to millionths with exact integers. The two bands are disjoint, so their squared norms ADD, rather than adding the norms and squaring. Dropping the parity restriction in the analytic remainder is conservative.

## 5. Executed exact comparisons and mixed positivity

| N | parity (sine indices) | rho_band upper | rho_tail upper | rho total upper | mu_N |
|---:|---|---:|---:|---:|---:|
|172|odd|0.303771|0.151364|0.455135|0.25|
|172|even|0.210782|0.140924|0.351706|0.25|
|242|odd|0.153093|0.144773|0.297866|0.6|
|242|even|0.135713|0.133840|0.269553|0.6|
|300|odd|0.108180|0.140713|0.248893|0.75|
|300|even|0.103096|0.129484|0.232580|0.75|

These decimals are exact rational outward upper budgets, not eigenvalues. Weighted Cauchy-Schwarz gives ||Q_N Tp||² <= rho sum d_i|z_i|² (separately per parity). For any eta>rho with eta<mu_N, Young's inequality gives

    W(p+q,p+q) >= sigma||p||²+(1-rho/eta)sum d_i|z_i|²+(mu_N-eta)||q||².

At N=242 use eta=3/10, since both upper budgets are strictly below 3/10. At N=300 use eta=1/4, since both are below 1/4. With the previously strict analytical tail margins, both imply strict inequalities for all nonzero p+q, p in E14, q in Q_N H0 intersect D_log:

    N=242: W(p+q,p+q)>10^-14||p||²+(3/10)||q||².
    N=300: W(p+q,p+q)>10^-14||p||²+(1/2)||q||².

Thus at cutoff 242 the mixed positive sector has codimension 228 within H0, compared with the former 14,986. At cutoff 300 it has codimension 286. This is not a positivity certificate on all of E242 or E300, nor is it an upper bound on the number of negative directions smaller than the 172 already obtained from Q172 alone. The 172 mixed test is inconclusive: the present upper budgets exceed 1/4, which does not establish failure of the true Schur complement or a negative direction. The actual full variance matrix is not evaluated.

## 6. Reproduction and next step

From the package directory run:

    python cr_alp_041_energy_conditioned_schur_audit.py --verify
    python cr_alp_042_support_aware_prime_tail_audit.py
    python cr_alp_043_exact_band_audit.py
    python cr_alp_043_check_result.py

The first two predecessor checks and the new interval sweep were executed successfully. The interval sweep is the substantial computation; the small final checker independently compares the saved rational budgets and tail margins, not the underlying coefficient construction. Optional numerical-only pilot: `python cr_alp_043_signed_band_diagnostic.py` (NumPy/SciPy). The exact sweep was restarted after a time limit and accelerated by interval rotations; the completed run covers every index in both sectors. No fresh Lean run or independent full mathematical formalization is claimed. All continuum identities and input matrix enclosures remain inherited analytical premises.

Next CR-ALP-044: replace the trace comparison at 172 by a tighter matrix-valued bound with a rational fitted pole subtraction and sharper oscillatory/parity tail estimates, or use the 242 mixed theorem as the working cutoff for a hierarchical certificate of its 228 missing directions. Do not infer positivity by summing separately positive sectors. Complete fixed collar, arbitrary-support extension, and RH remain open. No receiving acceptance ledger is promoted.
