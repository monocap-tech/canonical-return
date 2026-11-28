# CR-INT-T1D-AN — Local arithmetic comparator and adaptive Neumann bracketing

**Disposition:** analytical comparison with executed integer/rational audit, conditional on the original fixed-window Weil/Gamma identity. Not Lean, independent review, remote CI, low/mixed positivity or RH. AM's old global high gate and AL's mixed gate remain separate.

## 1. Preserve the source

On I=(-2,2), q0=Gamma-A on the smooth doubly Mellin-pole-null source. A is AM's original sum of 24 prime powers and 48 signed shifts. Use its certified shift floors and upper coefficients (am_prime.py), and its exact integer supersolution weights a_i on 16384 physical cells. Define positive weights h_i=floor((H*a_i^5)^(1/6)), H=10^8; the sixth-power inequalities are checked exactly. With h=0 outside I, the weighted Schur inequality gives |<Af,f>|<=int_I V|f|^2, V=Ah/h. On cell i the audit bounds V by the upward-rounded integer millionth of sum_s u_s max(h_(i+j_s),h_(i+j_s+1))/(10^9*h_i). This includes every translated position and all intercell interactions. Flattening the AM weight is a new proposal; its global max need not improve AM.

## 2. Direct-sum Neumann comparison

For any finite interval partition I=union I_j, let L_j be the Neumann Laplacian on I_j and v_j a certified upper bound on V there. For a>0 the compressed whole-line resolvent is dominated in quadratic-form order by the direct sum of Neumann resolvents (a^2+L_j)^-1. Proof: the variational supremum 2 Re<f,u>-int_R(|u'|^2+a^2|u|^2) over H^1(R) is bounded by the sum of the independent suprema over H^1(I_j), dropping exterior energy and matching constraints. The digamma series m_Gamma(ξ)=m_Gamma(0)+sum_(k>=0) (2/a_k) ξ²/(a_k²+ξ²), a_k=2k+1/2, has positive summands, so monotone summation gives Gamma>=direct sum m_Gamma(sqrt(L_j)). Consequently q0[f]>=sum_j <(m_Gamma(sqrt(L_j))-v_j) f_j,f_j>. This retains original physical L2 mass, not an identification of the local comparator with the original operator. The internal Neumann boundary defect is positive. For smooth f, restrictions belong to each local form domain.

If r_j is a positive integer with m_Gamma(pi*r_j/|I_j|)>v_j+epsilon, the local modes n>=r_j are positive with margin epsilon by monotonicity. Impose orthogonality to all local Neumann modes n=0,...,r_j-1. This excludes at most R=sum r_j dimensions, even on the pole-null source: any negative-definite subspace of larger dimension intersects the positive sector nontrivially. Reflection-symmetric partitions and paired r_j split R equally between parities. No additional subtraction for previous low modes or pole constraints is taken.

## 3. Improved exact Gamma lower floor

With x=|ξ|/2, t_k=k+1/4, T=129/4, write F_x(t)=x²/[t(t²+x²)]=1/t-t/(t²+x²). It is positive decreasing and convex: F''=2/t³-2t(t²-3x²)/(t²+x²)^3>=0 (if t²>=3x² compare denominators; otherwise the second term is nonnegative). Convex trapezoid comparison on unit intervals implies sum_(k>=32) F_x(t_k)>= integral_T^infty F_x(t)dt+F_x(T)/2. For x>T, integral=.5 log(1+x²/T²)>log(x/T). Also F_x(t)>=1/t-t/x². Thus a rigorous lower bound is m0+sum_(k<32)1/t_k-sum_(k<32)t_k/x²+log(x/T)+1/(2T)-T/(2x²). Use AM's rational lower m0 and pi bounds and rational lower log; x may be replaced by a positive lower value since the displayed bound increases with x. The audit checks each local floor strictly above v_j+1/20000.

## 4. Executed results

Equal 256 intervals: potential upper thresholds range 7.393152 to 8.537422, local exceptional counts 51 to 160, total 21820 (10910 per parity). Adaptive partition: a deterministic floating DP proposes intervals up to 512 cells on half the source, mirrored on the other half. It proposes integer ranks using exp(v); all endpoint coverage, max potential and final Gamma inequalities are checked with exact integers/Fractions independently of the floating proposal. Final 924 intervals, total 21574 (10787 per parity). Same improved Gamma floor with uniform 7.98 comparator: certified N=23379. Thus localization reduces this comparable bound by 1805; total improvement from AM's 23744 is 2170. Both local audits PASS. Conditional on the analytical imports, n_-(q0;D^0)<=21574, with q0>||f||²/20000 on the local-high smooth pole-null sector. No claim of a global Neumann high cutoff at 21574; the exceptional basis is now piecewise Neumann.

## 5. Structural next gate

The local comparator is not itself positive on all pole-null vectors: its 924 piecewise constant modes have negative diagonal values m_Gamma(0)-v_j, and their intersection with the two moment constraints has dimension at least 922 in the comparator form domain. This is NOT an original negative witness. Restoring the discarded terms is essential. For each positive shift s=log n, the exact weighted graph energy identity is

 <Af,f> = int_I V |f|² - sum_n kappa_n int_{I intersect (I-s)} h(x)h(x+s)|f(x)/h(x)-f(x+s)/h(x+s)|² dx.

Therefore the original form is the local diagonal comparator plus three nonnegative forms: the complete Gamma resolvent bracketing defect, the weighted arithmetic graph energy, and int (v-V)|f|². This exposes an exact possible positive correction mechanism rather than an uncontrolled subtraction. AO should test whether retaining the graph energy and Gamma matching corrections can control whole classes of local constant and low modes; do not claim positivity from the diagonal comparator alone. No mixed relative rerun at the new partition was executed. Complete fixed-window and global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.

Reproduce: python3 an_local.py (NumPy proposal, exact integer/Fraction checks). The AM scripts are included unchanged; audit.log is the executed output.