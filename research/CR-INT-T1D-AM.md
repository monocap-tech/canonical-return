# CR-INT-T1D-AM — Fixed-cell joint prime supersolution and reduced high cutoff

**Disposition:** new analytic weighted Schur comparison, with executed integer and rational audits. Conditional on the original fixed-window Weil/Gamma identities and AE's Neumann resolvent comparison. Not Lean, remote CI, independently reviewed, complete fixed-window positivity, all-support positivity, Suzuki–Connes transfer, Phase IV or RH.

## 1. Preserve the source

On I=(-2,2), keep all 24 original prime powers n<=54 and 48 signed shifts. A is the bounded self-adjoint positive-kernel auxiliary operator sum kappa_n P_I(T_log n+T_-log n)P_I*, kappa_n=Lambda(n)/sqrt(n). The signed Weil term is -<Af,f>. Positive-kernel does not mean positive semidefinite. Q previously proved ||A||<203/25 by a path-event supersolution. AL retained this norm at N0=27400, with conditional negative-index bound 27294 after its low/high mixed certificate.

## 2. A fixed-cell supersolution without point sampling

Partition I into M=16384 equal half-open cells of length delta=4/M (endpoints are immaterial). Let h be constant h_i>0 on each cell, and define h_i=0 for i outside 0,...,M-1. For each signed shift s=+/-log(n), certify the unique integer j=floor(s/delta) using atanh-series rational intervals; s/delta is strictly between j and j+1. On cell i, h(x+s) is bounded by max(h_(i+j),h_(i+j+1)). This includes zero extension and every possible translated position, not just a sample point. Let u_s be an integer upper bound on 10^9*kappa_n from the same log intervals and integer square root. Thus (Ah)(x)<=10^-9 sum_s u_s max(h_(i+j),h_(i+j+1)).

The audit proposes h by 80 deterministic integer iterations of that upper map, normalizing to max h<=10^8. The proposal uses NumPy int64 with verified bound 48*10^9*10^8<2^63 for nonnegative row sums, and unbounded Python integers for normalization. An independent final verifier using Python integers checks all 16384 inequalities

    50 sum_s u_s max(h_(i+j),h_(i+j+1)) < 399*10^9*h_i.

Executed PASS. min h=39499328, max h=100000000. Max upper ratio diagnostic 7.9796475323000635. For each shift the weighted inequality 2|f(x)f(x+s)| <= |f(x)|² h(x+s)/h(x)+|f(x+s)|² h(x)/h(x+s) applies on the overlap; symmetry of signed shifts gives |<Af,f>| <= integral (Ah/h)|f|². The finite maximum strict bound implies ||A||<399/50=7.98. This is a continuum bound without discretizing the actual quadratic form.

## 3. New high Neumann gate

AE established Gamma >= m_Gamma(sqrt(L_N)) in form sense. Its boundary correction is nonnegative. m_Gamma is increasing in absolute frequency by the positive digamma series. For n>=N, x=n*pi/8, t=129/4, the first 32 summands plus decreasing integral tail give

    m_Gamma(n*pi/4) > m_Gamma(0) + sum_(j=0)^31 x²/((j+1/4)((j+1/4)²+x²)) + log(x/t).

The standalone am_high.py audit uses Machin alternating pi bounds, harmonic bounds on Euler's constant as in R, atanh log bounds and exact fractions, rounding down. Subtract 399/50. The certified floors at N=23744,23800,24000,27400 are 53/10^6,2409/10^6,10777/10^6,143268/10^6, respectively. Hence q0[y]>||y||²/20000 on the pole-null high smooth sector n>=23744. The low Neumann projection has rank at most 23744 (11872 per parity), giving n_-(q0;D^0)<=23744 by the dimension argument. This is a reduction of 3550 from AL's conservative conditional bound; it does not assume the new cutoff can be combined with AL's V106. No conclusion about actual negative index.

## 4. Mixed rerun: incomplete, no promotion

An exploratory rerun of AL at N0=24000 and new arithmetic norm used its boundary threshold j>=8000, the same signed high Gram through 48000, and adjusted the exponential discrepancy check to sum exp(-4a_j)<1/20000 and w>10000 (so 2 sum exp/w<10^-8). Both high Gram computations completed, with leading squared upper diagnostics .001103446589499251 and .0011007545541739372 (old AL .001310404275673795 and .0013035506089984765). The two runs timed out before completion of the low structured factorization and final mixed sign tests. These high Gram diagnostics are not a relative beta certificate and no V106 plus lowered high-sector positivity is asserted. AL's earlier certified statement at N0=27400 remains distinct. The partial logs and scripts are retained explicitly as exploratory material.

## 5. Reproduction and next gate

Run `python3 am_prime.py` (NumPy for deterministic proposal; exact Python integer final verifier), then `python3 am_high.py` (standard library only). See prime.log and standalone_high.log for executed results. The exploratory am_run.py and dependencies are not part of the passed gate. No negative original witness. Remaining intermediate-sector sign and complete source transport OPEN.

**Next AN:** use the improved prime bound to complete the relative mixed rerun at a practical N0 such as 24000, or derive a local arithmetic comparator using the supersolution's pointwise Ah/h rather than its global maximum. Do not count an incomplete mixed rerun as a passed certificate; retain AL's original low certificate separately.