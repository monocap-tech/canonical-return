# CR-INT-T1D-AJ — Infinite positive Gamma boundary tail in the high inverse

**Disposition:** DERIVED / ANALYTIC; executed exact-rational interval audit PASS for 82 directions (AI comparison) and 90 directions (new certificate). Conditional on the original AE–AI form identities and Q's joint arithmetic norm. Not independently reviewed, Lean-certified, remote CI, whole fixed-window positivity, Suzuki–Connes completion, Phase IV, or RH.

## 1. Source and high comparator

Keep the physical interval I=(-2,2), both Mellin moment constraints, all 24 signed prime-power interactions, and the complete Gamma multiplier. Use AI's high cutoff N0=27400 and its 53 rational nondecreasing lower bands for D_h(n)=m_Gamma(nπ/4)−203/25. In each parity σ, retain the exact positive Gamma terms B_j for j≥K0=8000 in the high lower comparator. Thus on the pole-null high form space Y, q0[y]≥d[y]+b[y], with d[y]=<D_h y,y>, b[y]=Σ_{j≥K0} |<u_j,y>|². In normalized high cosine coordinates, a_j=2j+1/2, w_n=nπ/4, and u_j(n)=a_j sqrt(1−σ exp(−4a_j))/(a_j²+w_n²), up to the common parity sign, which can be absorbed in the rank-one factor. This is exactly AE's boundary correction, not a replacement by a periodic or truncated Gamma operator. B is interpreted as a positive form; no boundedness of the full infinite boundary operator is assumed.

For any ℓ² sequence z of boundary dual coordinates for which v=Σ u_j z_j exists in the weighted high dual space, Cauchy–Schwarz on the direct sum gives

    |<r,y>|² ≤ (||D_h^(-1/2)(r−v)||²+||z||²)(d[y]+b[y]).

Equivalently this is an upper bound for the inverse quadratic form of the strengthened comparator (with restriction to Y). It does not require construction of its full inverse. An unprojected representer can be used for this upper bound; AI's addition of a high pole row is retained because it pairs to zero on Y.

## 2. Exact cancellation of the retained boundary component

For low raw cosine coefficients c_m, set S_k=Σ c_m w_m^(2k), T_k=Σ c_m Fσ(w_m) w_m^(2k), as in AI. The raw high mixed leading coefficient is S0 times

    Fσ(w)/w² = 2Σ_{j≥1}(1−σ exp(−4a_j))/(a_j²+w²) + 2Σ κ_h sin(w log h)/w.

Choose z_j=S0 sqrt(2(1−σ exp(−4a_j)))/a_j for j≥K0. Then the normalized high vector v cancels exactly the j≥K0 portion of the first term (raw-to-normalized conversion 1/sqrt(2)). This series converges coordinatewise and in the weighted high dual norm, by comparison with the full decreasing integral majorant O(1/w); z is square summable. Its cost satisfies

    ||z||² ≤ |S0|² (10001/10000)/(2(K0−3/4)).

Here exp(−4a_j)<1/10000 and Σ_{j≥K0} a_j^(-2)≤1/[4(K0−3/4)] by the integral test. Only the residual Gamma leading amplitude remains:

    gσ(w)=2w Σ_{1≤j<K0}(1−σ exp(−4a_j))/(a_j²+w²).

It is positive and decreasing for w≥w_N0>2K0, because every summand is. The elementary decreasing integral test bounds it between atan((2K0+1/2)/w)−atan((5/2)/w)−10^(-8) and atan((2K0−3/2)/w)−atan((1/2)/w)+10^(-8). The exponential error uses (2/w)Σ_{j≥1}exp(−4a_j)<10^(-8). Both arctan arguments are less than one, so alternating rational series give directed enclosures.

## 3. Actual weighted mixed bound

For the first 12 bands (27400≤n<28000), evaluate the full square of gσ(w_n)+2Σ κ_h sin(w_n log h) with the exact interval complex recurrence of AI, enclosing gσ using its monotone endpoint bounds. For n≥28000 use the same integral and Abel majorants as AI; the new cross term is controlled by monotonicity of gσ(w_n)/(n² d(n)). All 24 phases are retained. For the direct-sum leading coefficient, the squared normalized bound is the resulting weighted square divided by two plus the displayed boundary dual-coordinate cost. The remaining S_k,T_k terms (k<24), the exact pole quotient of T0, and geometric remainder are unchanged, with the normalized coefficient bound 1/sqrt(2)<707107/1000000. Triangle inequality in the direct-sum space bounds the full mixed functional relative to the actual low energy, using the directed LDL energy-dual norms.

For AI's 41 dimensions per parity the 10^−220 implementation passes. Upper diagnostics for the new squared relative penalty: 0.8156170451162552 (even), 0.7953147597908526 (odd), compared with AI's 0.8881874986205489 and 0.8590726339711204. Exact assertions retain ratio<9/10. Boundary coordinate cost upper diagnostic is 0.00006251211051036 per |S0|²; residual weighted square upper diagnostics are 0.0013470383381303291 and 0.001323205162265041. These are upper bound diagnostics, not actual eigenvalues or precise improvement of the true inverse.

## 4. New 90-dimensional certificate

Increase AI's raw parity dimension from 42 to 46 (45 exact pole-null vectors per parity; raw modes p,p+2,...,p+90). Retain the exact same physical metric, complete low matrix and N0. The direct generator interval LDL is run on a 10^−320 grid, normalized Euler–Maclaurin parameters P=110, K=500, R=150, and correspondingly strengthened Machin, atanh, exponential and trigonometric series. Both independent parity runs terminate PASS. All 90 pivots are positive. Reciprocal physical trace lower diagnostics: 1.400724620926121e−92 (even), 7.553924829198637e−90 (odd); exact assertions give q0[f]>10^−93||f||² on V90. Relative mixed bound upper diagnostics: 0.948285828770256 and 0.9106894775644911. Exact assertions in both runs give

    |q0(f,y)|² ≤ (97/100) q0[f] (d[y]+b[y]).

Since (197/200)²>97/100, the original complete form obeys

    q0[f+y] ≥ (3/200)(q0[f]+d[y]+b[y])

on V90⊕Y. In particular it is strictly positive there for nonzero vectors, conditional on Q and AE–AI. The same conservative codimension counting as AI gives n_−(q0;D⁰)≤27400−90=27310. (The separate pole-null constraints are not subtracted from this conservative count.) This does not constrain the actual count more tightly than that upper bound.

The attempted 49-dimensional even low block passed the interval LDL but its mixed upper diagnostic 1.0790542794141618 did not pass the positivity threshold; no corresponding odd result or negative original direction is certified. The baseline 58-raw-mode attempts failed interval pivot resolution. Do not confuse an upper majorant exceeding one, or an interval enclosing zero, with a negative spectral witness.

## 5. Reproduction and next gate

Bundle: baseline `python3 aj_boundary.py` with the three AI dependency modules; high-precision `cd high_precision; python3 run90.py 0; python3 run90.py 1` (separate parity runs). Logs `audit_output.txt`, `audit90_even.txt`, `audit90_odd.txt` record executed tests. `explore50.txt` and `explore58.txt` document unsuccessful exploratory gates. Integer intervals and Fraction assertions are exact; floats only print diagnostics. The written analytical source identities, variational inequality, infinite series and tail estimates are not mechanically verified by the scripts.

**Next AK:** replace the leading-only boundary dual ansatz and separate moment Minkowski bounds with a joint energy-dual Gram certificate, using additional boundary coordinates and signed finite-band correlations. Test whether a useful larger block can pass without moving the high cutoff; alternatively tighten the original joint arithmetic norm. The remaining intermediate sector and mixed interactions, original pole directions, support-independent positivity, completed Suzuki–Connes transfer, Phase IV and RH remain open.