# MASTER-CROSS-004 — prime-channel globalization obstruction

2026-11-26 PST. INVESTIGATE / ADOPT 0 / NONPROMOTING. Continues custody commit b26e48c530c1cef50f4b0626c0dc391f289ece3d, blob e100790dfeae414c458e9dc73baa7c75967be83a. New human analytic deductions, no Lean compile, no numerical positivity experiment or RH proof.

## 1. Reopened source contracts

Canonical-return d14933facc14e15cb9be07f09ce36f5eec9b83be:
- research/CR-ALP-024_FIXED_COLLAR_WEIL_LOG_SPECTRAL_FINITE_NEGATIVE_INDEX.md d9dc473a269ac5f2ed7be6e403699e8cb5d094f0: signed original W=Gamma+pole-prime, pole-null reduction, log(1+|r|)-8 <= g(r) <= log(1+|r|)+3.
- research/CR-ALP-042_SUPPORT_AWARE_PRIME_POWER_NORM_AND_172_MODE_TAIL.md 395c949e6fda0d726a757cc8bfa281370273eb89: seven-term positive-weight self-adjoint translation operator, fixed collar bound 21/8, not a uniform support bound.
- research/CR-ALP-029_ORIGINAL_SIGNED_SIX_MODE_GAMMA_CERTIFICATE.md e4c08bc0f192a7874622928b6ee03f6545b753d2: Gamma form = g(0) mass + sum_{m>=0} integral_0^infty exp(-beta_m s)||f-tau_s f||² ds, beta_m=2m+1/2, inherited classical digamma identity. Polarization gives mixed resolvent form.

External HMarcusWH/zeta-23-lean pinned 7597acdba0c7d0a81b8c3a5f229598e3f9b74ac9; Step123NativeFieldRecord bff5d4d2011a0d0da973605033acbcfc59f67a4b, Step123OutputData 34e33b588f1b613b97d3dc436f71b5c46e71d659, Step123CompleteOutput 65086c6ca38767c6ef524b7f66e8b58658976460 as CROSS-003. Round6 README 6debff2b5a382e30508aa996374138688ce5389b directly reopened: native finite signed path identities do not supply the full growing-cutoff upper; independent ACTUAL_UPPER and PAPER_ROUTE_VALIDATION remain open. Round4 closure register 22928216dafe6574fe302c896bfeebe134d0894b includes C2_ACTUAL_LONG_TERM_LOWER and consumer riemannHypothesis_of_actualCorrelationLongTerm_eventual_lower. Exact C2 consumer source not located by code search in this turn; do not assert equivalence of its literal source to the new scalar formula below.

## 2. Prime-only spectral-edge invariance (new analytic lemma)

For any A>0, I=(-A,A), define B_A=sum_{log n<=2A} w_n P_A(tau_{log n}+tau_{-log n})P_A, w_n=Lambda(n)/sqrt(n)>=0, tau_t f(x)=f(x-t), P_A zero-extended physical projection. B_A bounded self-adjoint and positivity preserving, not necessarily PSD. Let b(A)=sup_{||f||=1} <f,B_A f>. Since |<f,B_A f>|<=<|f|,B_A |f|>, b(A)=||B_A||. For either centered physical parity sigma=±1, let S_{A,sigma}^0 be compact smooth tests in I with that parity and both F_±=integral exp(±x/2)f(x)dx=0.

Claim: sup_{f in S_{A,sigma}^0,||f||=1} <f,B_A f> = b(A) for BOTH sigma. This holds for every A, even if prime inventory empty.

Proof:
1. Approximate the upper edge by nonnegative even h of unit mass (L2 norm one). For arbitrary f, put h(x)=sqrt((|f(x)|²+|f(-x)|²)/2). For each paired translation form, pointwise two-vector Cauchy-Schwarz and reflection symmetry give <h,B_A h>>=<f,B_A f>; same physical mass. Approximate h in L2 by even nonnegative compact smooth h and normalize. B_A bounded, so Rayleigh quotient approximates the edge.
2. For the finite shifts ell_n=log n, choose t_k->infty with exp(i t_k ell_n)->1 simultaneously. Compactness of the finite torus applied to integer-time orbit, taking sufficiently widely separated pairs in a convergent subsequence, gives positive differences tending to infinity and phases tending to identity. No prime log independence is needed.
3. Define v_{k,+}=sqrt(2)h cos(t_k x), v_{k,-}=sqrt(2)h sin(t_k x). Norms tend one by Riemann-Lebesgue. For every shift ell, their self-correlation equals cos(t_k ell) times h autocorrelation plus/minus an oscillatory integral at frequency 2t_k, tending to zero. Thus <v,B_A v> tends <h,B_A h> for both.
4. F_+(v_k)->0 by Riemann-Lebesgue. For same parity F_-(v_k)=sigma F_+(v_k). Choose once a real compact smooth parity-corrector c_sigma with F_+(c_sigma)!=0 (even positive bump, odd x times even nonnegative bump). Subtract c_sigma F_+(v_k)/F_+(c_sigma). Both poles vanish exactly, parity and support unchanged, correction tends zero in L2. Normalize. Bounded B_A preserves limiting Rayleigh quotient. Since initial h approached b(A), result follows. No signed Weil positivity claim.

Growth: normalized indicator of I gives exact b(A)>=2 sum_{log n<=2A} w_n(1-log n/(2A))>=sum_{n<=exp(A)} w_n. Diverges as A->infty by divergence of sum_p 1/p (w_p>=1/p). Thus no aperture-independent mass-only upper on positive prime channel even on either smooth pole-null parity. Pole projection and symmetry cannot improve its upper spectral edge at any aperture. For fixed A the oscillatory corrected tests also have Gamma log energy tending to infinity (low-frequency mass tends zero by Fourier shift and correction norm->0). Hence not negative Weil witnesses; cannot infer any sign of complete W. The result rules out uniform prime-only norm strategy, not ALP fixed collar or Gamma-relative comparisons.

## 3. Exact continuous main-term subtraction

Define C_A on physical interval with symmetric kernel exp(|x-y|/2), K_A with kernel exp(-|x-y|/2). Then C_A=(exp((x-y)/2)+exp(-(x-y)/2)) kernel minus K_A. The first kernel is rank at most two and vanishes in quadratic form on the two-pole-null carrier (and compressed action vanishes). Thus <f,C_A f>=-<f,K_A f>. K_A positive semidefinite by Fourier transform 1/(r²+1/4) of exp(-|t|/2), and ||K_A||<=4 by convolution L1 norm. Put E_A=B_A-C_A. Exact identity on pole-null carrier:
W_A(f,f)=G(f,f)+<f,K_A f>-<f,E_A f>.
The signed arithmetic defect is discrete prime measure minus exp(|t|/2)dt, not a positive kernel. It retains unbounded upper spectral edge on either parity because C_A compressed norm<=4. So even continuous main-term subtraction cannot supply uniform mass-only bound on E_A. Main term here is a chosen exact subtraction, not an approximation assumption.

By global Gamma estimates and ||K_A||<=4, the target E_A[f]<=E_log[f]+C||f||² uniformly on one parity is equivalent (up to additive constants) to corresponding uniform W_A>=-C'||f||². Likewise B_A[f]<=E_log[f]+C||f||² is equivalent; Gamma-log difference bounded. Coefficient 1 is the natural exact comparison; no assertion that it is optimal or a bound exists. Under CROSS-003 obstruction supply, such uniform relative inequality implies RH, so cannot be counted as already established progress toward it.

## 4. Narrower source-specific scalar receiver

External original trial is a fixed real even compact two-pole-null kernel g proportional to k=Pq, P=d²-1/4, normalized ||g||2=1, support [-b,b] with b=1/4096. Use the original fixed kernel, not an arbitrary chosen bump. For t>2b let u_{t,±}=(tau_{-t/2}g ± tau_{t/2}g)/sqrt(2); disjoint unit lobes, same two pole constraints and centered parity. These are the centered original RHRC trials at separation t=2a, L=t+2(1+b). Define autocorrelation h(s)=int g(x)g(x-s)dx, real even supported [-2b,2b]. Put
P_g(t)=sum_{n>=2} Lambda(n)/sqrt(n) h(t-log n).
At each t finite, includes every prime power in the support, no absolute values. For t>2b the opposite-shift terms h(t+log n) vanish. Translation invariance and polarization give the exact real diagonal:
W(u_{t,±},u_{t,±})=W(g,g) ± [Gamma_cross(t)-P_g(t)].
The cross Gamma term is controlled from ALP-029: polarized J_beta=2/beta identity minus convolution by exp(-beta|x-y|). For disjoint lobes the identity term is zero. The beta_0=1/2 exponential kernel contribution factors into the two elementary pole moments and vanishes exactly. Therefore, for d=t-2b>0,
|Gamma_cross(t)| <= ||g||_1² exp(-5d/2)/(1-exp(-2d)).
No truncation error; absolutely convergent resolvent expansion at separated supports.

Continuous comparison cancels exactly in this scalar source: integral_0^infty exp(s/2)h(t-s)ds=exp(t/2)int h(v)exp(-v/2)dv=0 for t>2b, since latter factorizes into F_-(g)F_+(g). Hence P_g(t) is an exact signed prime discrepancy against the chosen continuum main term.

Under an off-line zero, original_depth and cofinal separate phase supply give, for fixed 0<beta<zero displacement, sequences t_+,t_- tending to infinity with P_g(t_+) >= 4 exp(beta (t_++2(1+b)))+W(g,g)+Gamma_cross(t_+) and P_g(t_-) <= -4 exp(beta (t_-+2(1+b)))-W(g,g)+Gamma_cross(t_-). These follow by rearranging exact trial identities. Consequently either eventual upper bound P_g(t)<=C, or eventual lower bound P_g(t)>=-C, for all sufficiently large real t suffices to exclude off-line zeros. More generally one-sided subexponential bound, for every eps>0 with constants fixed before t, also suffices by choosing eps<beta. Merely cofinal t upper bound is not sufficient from these separate cofinal witness sequences. No bound of this strength established. This is likely close to external C2 long-term route; exact identity/redundancy not yet checked. The reduction is conditional on original kernel identification and inherited analytic Gamma expansion, no Lean compilation.

## 5. Governance and next

T3 strengthened: generic pole-null, parity and continuous-main-term prime norm routes rejected for uniform aperture control. T1 analytic fixed-kernel scalar reduction, source-specific; formal attachment pending. T4 compare with external C2 actualCorrelationLongTerm consumer and RECON source arithmetic before declaring distinct. T5 retain fixed-collar ALP and route signed prime discrepancy to RECON / RETURN comparative interface. Round6 signed survivor charges are different observables; no transport asserted merely because both involve arithmetic cancellations. Next CROSS-005: locate exact external C2 source, match scalar long term including sign/normalization/quantifiers, inspect existing RECON weighted sampling and signed prime-discrepancy transport, test one-sided subexponential supplier. No arbitrary-support positivity or RH proof. No workflow explicitly launched. Commit note to existing isolated custody branch; do not merge or promote.
