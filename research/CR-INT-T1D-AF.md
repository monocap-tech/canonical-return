# CR-INT-T1D-AF — Certified near-null parity Ritz directions and a Gamma truncation obstruction

**Disposition:** DERIVED / exact rational interval AUDIT PASS, conditional on the AE source identity, boundary-resolvent representation, and form-core transfer. No Lean, remote CI, independent review, unconditioned positivity, completed Suzuki–Connes transfer, Phase IV, or RH.

## 1. Source, parity and physical metric

Retain AE's original signed pole-null Weil pairing on I=(-2,2), with the two Mellin constraints and unitary Fourier transform. On the form completion it is q₀=Γ−A, with Γ=D+Σ B_j and all 24 true prime powers included. Write y=x+2, ω_n=nπ/4, p_n(y)=cos(ω_n y). These are *unnormalized* Neumann functions: physical mass is 4 at n=0 and 2 otherwise. Reflection parity is (-1)^n.

For σ=±1 let n_i=2i+p (p=0 even, p=1 odd), i=0,...,7. The same-parity moment ratios against either exponential are exactly

    r_i=(1/4+ω_{n_0}²)/(1/4+ω_{n_i}²), i=1,...,7.

Use f_i=p_{n_i}−r_i p_{n_0}. Their two Mellin moments vanish exactly, including in the odd sector. Physical Gram G_ij=2 δ_ij + d₀ r_i r_j, with d₀=4 (even) or 2 (odd). These functions belong to the logarithmic form completion, not literally C_c^∞. AE's form-core density with exact moment correction supplies smooth pole-null approximants. Do not confuse finite positivity with a theorem on the full smooth core.

## 2. Complete original matrix and certified arithmetic

The Gamma matrix before congruence is diagonal d_n mΓ(ω_n) plus the same-parity entries

    B_j[m,n]=2 a_j² (1−σ exp(−4a_j))/[(a_j²+ω_m²)(a_j²+ω_n²)], a_j=2j+1/2.

The j=0 quadratic term vanishes exactly after the pole-null congruence, not off that space. For s=log h, κ=Λ(h)/sqrt(h), the sum of the two compressed shift overlaps in the unnormalized cosine basis for distinct same-parity indices is

    [(sin(ω_m s)−sin(ω_n s))/(ω_n−ω_m)] − [(sin(ω_m s)+sin(ω_n s))/(ω_n+ω_m)].

For n=m>0 it is (4−s)cos(ω_n s)−sin(ω_n s)/ω_n, and for n=m=0 it is 2(4−s). Multiply each by κ and sum over all 24 prime powers ≤54. Subtract this arithmetic matrix from Γ. These formulas follow by product-to-sum and the fact that (n±m)π is an even multiple of π. Opposite parity cross entries vanish. Congruence by f_i yields the physical q matrix Q, retaining its true signs and G.

The standard-library audit computes all intervals on a directed rational 10^−55 grid. No floating-point comparisons are used in the assertions. Machin's arctan series gives π; the atanh series supplies logs; integer roots bound sqrt(h); range-reduced Taylor recurrences with explicit remainders give sine/cosine. A positive Taylor series with a geometric remainder bounds exp(−u) by reciprocal intervals.

For Gamma sums set t=64+1/4, x=ω/2. The digamma expansion is evaluated by the first 64 terms and 22 geometric powers of x²/t². For the infinite tail use Hurwitz sums Z_s(t)=Σ_{j≥0}(t+j)^−s. Euler–Maclaurin at order 2P (P=16) gives

    Z_s(t)=t^(1−s)/(s−1)+t^(−s)/2+Σ_{k=1}^P B_{2k}(s)_{2k−1} t^(−s−2k+1)/(2k)! + E,
    |E| ≤ |B_{2P}| (s)_{2P−1} t^(−s−2P+1)/(2P)!.

This remainder follows from the bound on the periodic Bernoulli polynomial by |B_{2P}| (its absolutely convergent Fourier series). The harmonic Euler–Maclaurin formula bounds Euler's constant similarly at N=128. For mΓ the alternating geometric tail remainder is at most x^(2R+2) Z_(2R+3)(t), R=22.

For the B series, sum j=1,...,63 explicitly. With X=ω_m²/4, Y=ω_n²/4, the nonexponential tail equals

    (1/2) Σ_{r≥0} (−1)^r Σ_{k=0}^r X^k Y^(r−k) Z_(2r+2)(t).

For ρ=max(X,Y)/t²<1 the absolute error after r=R is bounded by

    (Z_2(t)/2) ρ^(R+1) [(R+2)−(R+1)ρ]/(1−ρ)².

The omitted exponential part is bounded by exp(−8·64−2) Z_2(t)/2. These bounds control the entire infinite boundary correction, not only a finite positive truncation. The finite matrix entries and pole-null congruence use directed intervals. Interval symmetric elimination of Q−10^−25 G certifies every pivot strictly positive in both parities. Therefore on the direct sum of these 14 trial dimensions q₀[f]>10^−25||f||² for nonzero f. The finite trial space is in the completion. No positivity of its infinite complement follows.

## 3. Near-null witnesses and physical spectral upper bound

The auditor computes midpoint rational elimination on Q to propose a vector v with last coefficient 1 and Q's first six rows approximately annihilated. It rounds the coefficients to the 10^−20 rational grid, then independently encloses vᵀQv and vᵀGv by directed intervals. This construction does not use midpoint diagonalization as a proof. The actual Rayleigh quotients satisfy the exact rational comparisons

    even: 87665/10^26 < q₀[f]/||f||² < 87670/10^26,
    odd:  73545/10^25 < q₀[f]/||f||² < 73557/10^25.

Thus the full smooth pole-null Rayleigh infimum is <10^−21 by the form-core density of AE, even if it were negative. This is an upper bound only, not a proof of nonnegativity, or an exact smallest eigenvalue. The result exposes a real small-energy direction rather than a machine-precision artifact. In particular no uniform lower bound ≥10^−21 is possible on the original fixed-window smooth pole-null form.

## 4. A rigorous obstruction to simple Gamma truncation

For the two same trial vectors write coefficients d_n of the unnormalized cosines, with d_{n_0}=−Σ r_i v_i. For j=32 evaluate the exact boundary energy

    B_j[f]=2a_j² (1−σ exp(−4a_j)) [Σ_n d_n/(a_j²+ω_n²)]².

The rational audit verifies B_32[f]/||f||² > 27/10^22 (even) and >16/10^21 (odd). Each is strictly larger than the full original q₀ Rayleigh quotient above. Since all B_j are nonnegative, D+Σ_{j<32} B_j−A is strictly negative on these same strictly positive original trials. The same applies to any K≤32; no negative direction of q₀ follows. The attempted stronger claim for K=10^6 did not pass the inequality check and is not asserted.

Consequently AE's sufficient Schur matrix using K≤32 cannot certify positivity globally, irrespective of low-block dimension or mixed-tail accuracy. This is a structural failure of that particular lower approximation: even before the high block is eliminated the truncated form has a negative constrained trial. Enlarging K without accurately restoring the analytic tail is not a sign certificate. Coarse absolute norm error allowances are particularly unsuitable near these tiny Rayleigh values.

## 5. Audit and next gate

Run `python3 canonical_return_t1d_af_neumann_ritz.py`. The attached output is the local run, with exact comparisons in code and displayed approximate diagnostics. The derivation imports the AE source/form identities and standard Euler–Maclaurin estimates; it is not kernel checked. No original negative direction was found. The full mixed low/high Schur matrix and its sign have not been certified. AE's conditional negative-index bound remains 32768. The only new positive certificate is the 14-dimensional pole-null compression.

**Next AG:** design a cancellation-resistant lower certificate that retains the *entire* Gamma boundary series analytically, or use a positive spectral-sampling comparison for the genuine small-eigenvalue sector with a certified remainder. A useful approach is to precondition by exact low Ritz Cholesky/LDL factors and certify *relative* mixed-tail errors in that energy metric, rather than subtracting a global scalar error. Establish whether the full low/high interaction satisfies the corresponding relative Schur bound. The finite certificate cannot be promoted to global positivity. Keep research draft; full fixed-window pole-null positivity, both physical pole directions, Suzuki–Connes transport, Phase IV and RH remain open.