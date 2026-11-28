# CR-INT-T1D-BM — Exact signed-operator tail moment cancellation

**Status:** local directed-interval and exact integer residual audits PASS in both parities. Construct constraint-isomorphic lifts with four exactly cancelled operator-action moments. Squared infinite action tail O(N^-9) in unnormalized original BL lift coordinates. Conditional negative-index bound remains <=240 (121 even,119 odd); BK complementary coercivity margins 17/10000 and 53/10000 retained. No projected finite-head residual or full Schur lower-bound sign evaluated. Conditional on inherited original fixed-window (-2,2), double Mellin-pole-null form identities, AF Neumann representation, AM norm bound, BK complement and form/core transfer, BL constraints. No Lean, remote CI, independent mathematical review, or frequency scan.

## 1. Exact signed tail identity

Use the original unnormalized Neumann modes p_n=cos(n*pi*(x+2)/4) on I=(-2,2), w_n=n*pi/4, sigma=(-1)^parity, and L=Gamma-A. For finite f=sum_n c_n p_n with one reflection parity, let x=w_m exceed every supported w_n. Put a_k=2k+1/2, J(w)=w^2 sum_k (a_k^2+w^2)^-1, kappa_q=log(p)/sqrt(q) for the 24 prime powers q=p^j<=54, and d(w)=J(w)+w sum_q kappa_q sin(w log q). The complete nonexponential Gamma and arithmetic off-diagonal formulas from BL yield exactly

    <p_m,Lf> = 2 d(x) sum_n c_n/(x^2-w_n^2)
               - 2 sum_n c_n d(w_n)/(x^2-w_n^2)
               - 2 sigma sum_{k>=1} a_k^2 exp(-4a_k)/(x^2+a_k^2) sum_n c_n/(a_k^2+w_n^2).

The k=0 exponential term vanishes on the pole-null lift. For opposite parity the coefficient is zero. The diagonal Gamma term is zero for exterior m. Define M_l=sum_n c_n w_n^(2l) and H_l=sum_n c_n d(w_n) w_n^(2l)+sigma*(-1)^l sum_{k>=1} a_k^(2l+2) exp(-4a_k) sum_n c_n/(a_k^2+w_n^2). These are the exact coefficients of the first terms of the two rational denominators after combining the signed operator. Sequentially cancelling M_0,H_0,M_1,H_1 yields pointwise pairing orders x^-2,x^-3,x^-4,x^-5 instead of generic x^-1. M_0=M_1=0 enforces endpoint values and second derivatives zero; first and third derivatives already vanish for finite Neumann sums. Exponential expansion uses finite algebraic remainder, not an unjustified convergent expansion in a_k/x for unbounded k.

## 2. Constraint-preserving four-moment correction

BL has V (d=121/119) and V_n (four extra pole-null Neumann modes), exact U=V_n-V C, C=(F*V)^-1 F*V_n. BL stores rational C0 and certified ||C-C0||_F<=e_C. Normalize moment rows by 1,1/W,1/W^2,1/W^3 respectively for M_0,H_0,M_1,H_1, W=4096>all retained w_n. Let M=[M_o,M_n] denote their evaluations on [V,V_n], K=M_n-M_o C. Enclose each row via BL's directed scalar special-function evaluator (including all 24 prime powers), k=1..39 exponential terms with explicit geometric tail from k=40. Final moment grid 10^-12; exact Frobenius error <=23/10^12 in both parities. With center K0 and rational inverse proposal Z, certify

    ||Z K-I|| <= ||Z K0-I||_F + ||Z M_o0||_F e_C
                + ||Z||_F e_M(1+||C0||_F+e_C) = rho < 1/2.

All accepted products and Frobenius comparisons are exact integers/rationals; numerical inverse only proposes Z. Thus exact K invertible and X=K^-1 M_o gives W_lift=V-U X with F*W_lift=F*V and every one of the four moments identically zero. This is not the four-mode energy-minimizing corrector from BL: its four degrees of freedom are spent on moment cancellation. The full 125/123-dimensional finite extension remains positive by BL; this does not establish the sign of the full Schur complement. Executed results:

| parity | rho (diagnostic) | upper ||X||_F | upper ||c(.)||_(2->1) |
|---|---:|---:|---:|
| even | .34073636447000005 | 7135 | 114000 |
| odd | .10179119906490508 | 953 | 17000 |

The last bound includes the pole-null reference mode coefficient. It follows from ||X||_F <= (||Z M_o0||_F+||Z||_F e_M)/(1-rho) and triangle/Cauchy-Schwarz for coefficient vector [(I+C X)z;-Xz] and its reference coefficient. Exact fractional values in NPZ. Bounds large, especially even. No claim of optimal correction modes or useful energy-normalized conditioning.

## 3. Certified complete infinite action tail

The rational sum bound sum_q kappa_q<16 and the decreasing integral comparison J(x)<=1+pi*x/4 give |d(x)|<=18x for x>=1 and |d(w_n)|<=18W. For x>=2W, x^2-w_n^2>=3x^2/4. The four cancellations leave remainder w_n^4/[x^4(x^2-w_n^2)] and a_k^4/[x^4(x^2+a_k^2)] respectively. The latter is bounded using sum_{k>=1} a_k^4 exp(-4a_k)<1, checked by a geometric bound from k=1. Hence

    |<p_m,L W_lift z>| <= 74 W^4 A_sigma ||z||_2 / w_m^5,

where A_sigma is the exact stored coefficient 2->1 bound. The nonzero Neumann modes have mass 2; summing over all integer m>N (not merely same parity) and pi>3 yields

    (Q_N L W_lift)*(Q_N L W_lift) <= tau_N I,
    tau_N=(74 W^4 A_sigma)^2 (4/3)^10/(18 N^9),   3N/4>=2W.

Q_N is the physical orthogonal Neumann high-index projection. Exact scalar acceptance checks:

| N | even tau_N upper | odd tau_N upper |
|---:|---:|---:|
| 262144 | 1e-6 | 3e-8 |
| 1048576 | 4e-12 | 8e-14 |
| 16777216 | 6e-23 | 2e-24 |

These are upper bounds on squared L2 operator-action tails per unit original lift-coordinate vector, not the residual in energy-normalized coordinates, not a bound on the high-frequency projection of P_C L W_lift, and not a full Schur correction. No enumeration to these cutoffs is needed for the tail bounds. A large finite-head calculation would still be needed using this crude estimate.

## 4. Weighted inverse and projection-safe head-tail separation

Let P_C be BK's physical complement projection within pole-null space. On C, q>=mu||f||^2 by BK; also q>=sum_n (m_Gamma(w_n)-lambda)||f_n||^2 by the positive Neumann Gamma boundary and AM lambda=399/50. The BL scalar evaluator verifies m_Gamma(0)>-a, a=43/8; the resolvent representation makes m_Gamma nondecreasing for w>=0. Set theta=mu/[2(mu+a+lambda)]. Convex combination gives on C the lower diagonal form

    q[f]>=sum_n d_n ||f_n||^2,
    d_n=mu/2+theta*(m_Gamma(w_n)+a)>0.

For any physical y define its complementary dual seminorm ||y||_{C,*}^2=<P_C y,L_C^-1 P_C y>. By the form variational inverse bound, ||y||_{C,*}^2<=sum_n ||y_n||^2/d_n, and also <=mu^-1||P_C y||^2. With y_N=P_{<=N} L W_lift z, t_N=Q_N L W_lift z, triangle gives

    ||L W_lift z||_{C,*} <= mu^-1/2 ||P_C y_N|| + sqrt(tau_N/d_N)||z||.

Use min(mu^-1,d_N^-1) for the tail if desired. This separation is valid even when P_C and Q_N do not commute. The projected head norm can be calculated via BK's exact Gram and pole projection: ||P_C y_N||^2=||y_N||^2-(F*y_N)*G^-1(F*y_N)-|<g,y_N>|^2/||g||^2, with F columns pole-null and g the parity Mellin representer. No head evaluation here. Exact scalar audit verifies d_N>mu for N>=262144 by monotonicity; at N=1048576 d_N>.0019419 even and >.0060533 odd (about 12.4% smaller tail inverse penalty than scalar). Near-null energy resolution still requires stronger constants/relative metrics and/or an efficient finite-head evaluator. In particular, fast asymptotic decay alone does not settle the 240-dimensional sign problem.

## 5. Reproduction and next BN

All 45 entries of BL's inherited SHA-256 manifest verified. BM files: bm_moments.py, bm_residual.py, propose.py (heuristic exploration only), actual integer moment enclosures and inverse proposals, two moment logs, scalar residual log, inherited digest log. Run `OPENBLAS_NUM_THREADS=4 python bm_moments.py 0` and `... 1`, then `python bm_residual.py` from the BM directory with sibling BL input directory. Compact archive includes the directly consumed BL evaluator and four input NPZ files with their original hashes; full upstream certification requires BL and its inherited dependencies. No independent formal verification. The arbitrary 4096 normalization and use of the existing four BL extra modes were chosen for this experiment, not optimized. The proof is conditional on inherited analytical and numerical certificate chain.

**Next BN:** compute the actual projected finite-head residual in original signed-energy-preconditioned coordinates, preferably using rectangular mixed-action formulas and weighted/structured tail estimates. Optimize moment-null corrections jointly with the finite energy-minimizing response, rather than forcing four poorly conditioned directions to cancel four moments. Conditional index remains <=240. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.