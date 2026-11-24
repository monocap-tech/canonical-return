# CR-INT-T1D-AG — Relative-energy mixed-tail certification by the complete divided-difference symbol

**Disposition:** DERIVED / ANALYTIC with executed exact-rational interval audit PASS, conditional on AE/AF's original source/form identities and T1D-Q's joint prime norm. Not independent review, Lean, remote CI, whole low-block positivity, completed Suzuki–Connes transfer, Phase IV, or RH.

## 1. Frozen source and the actual mixed block

Use the complete original signed pole-null Weil form q₀=Γ−A on I=(-2,2), retaining the physical L² metric, both Mellin moments and all 24 prime powers. Let V=V₊⊕V₋ be AF's 14-dimensional pole-null form-completion space, with seven trial vectors in each parity. In one parity σ=(-1)^p use the unnormalized modes p_n=cos(ω_n(x+2)), ω_n=nπ/4, n=p,p+2,...,p+14, with physical masses 4 at n=0 and 2 otherwise. The AF basis is f_i=p_{p+2i}−r_i p_p, r_i=(1/4+ω_p²)/(1/4+ω_{p+2i}²), i=1,...,7. Let H be its actual full signed q₀ Gram, and G the physical Gram. AF's interval certificate H>10^(-25)G is recomputed by the AG audit.

Put N₀=32768 and let Q_h project onto the physical orthonormal Neumann modes e_n=p_n/sqrt(2), n≥N₀. Define Y to be the intersection of this high sector with the original pole-null logarithmic form domain. AE/Q give D_h=diag(mΓ(ω_n)−203/25)>1/8 and q₀[y]≥<D_h y,y> for y∈Y. The mixed representer Xf=Q_h K f is defined by the complete pole-free Fourier multiplier K (compressed to I); its pairing with y∈Y is the actual original Weil pairing. For fixed f∈V its zero extension is piecewise smooth with Fourier decay O(1/|ξ|), so the logarithmically growing multiplier image is in L². Thus these mixed coefficients represent q₀(f,y) on the actual form domain, without assuming a stronger domain-core theorem for high Neumann projections. Pole terms vanish on V⊕Y. No whole-operator L² boundedness is assumed.

## 2. Complete scalar divided difference

For each parity define the real scalar function

    Fσ(w)=bσ(w)+a(w),
    bσ(w)=2 Σ_{j≥1} (1−σ exp(−4 a_j)) w²/(a_j²+w²), a_j=2j+1/2,
    a(w)=2w Σ_{h≤54} κ_h sin(w log h), κ_h=Λ(h)/sqrt(h).

The j=0 boundary term is omitted here only because the low argument is pole-null, so this term's pairing with every high mode is identically zero. For distinct same-parity raw modes, the diagonal Γ term contributes zero, while the boundary sum over j≥1 is (bσ(w_n)−bσ(w_m))/(w_n²−w_m²) by elementary partial fractions. AF's exact compressed signed-shift formula gives the prime operator entry (a(w_m)−a(w_n))/(w_n²−w_m²). Hence the mixed pairing, after the exact low pole-null constraint, is

    q₀(f,p_n)=Σ_m d_m [Fσ(w_n)−Fσ(w_m)]/(w_n²−w_m²),
    f=Σ_m d_m p_m, n≥N₀ of matching parity.

This uses the *entire* Gamma boundary series and the original signed prime form. Opposite parity mixed coefficients are zero. Raw-mode identities here use the pole-free Fourier expression and the omitted j=0 correction; only their constrained pairing is the original q₀. The formula is not a replacement by a positive arithmetic operator.

For w≥0, bσ(w)≤(1+exp(−10))πw/2<8w/5: since a_j≥2j, use the decreasing integral test on Σ_{j≥1} w²/(4j²+w²). The original sum 2Σκ_h<25 is verified in AG, so |Fσ(w)|≤30w. On the low indices w_m<12; in particular |Fσ(w_m)|<360. The code evaluates bσ(w_m)=w_m² bc(w_m,0,p) by AF's rigorously enclosed full boundary series, with bc being the raw boundary Gram formula evaluated at the auxiliary argument zero (also valid as an algebraic scalar in odd parity).

## 3. Relative dual norms and convergent mixed-tail expansion

Write S_k(f)=Σ_m d_m w_m^(2k) and T_k(f)=Σ_m d_m Fσ(w_m) w_m^(2k). The exact geometric expansion is

    q₀(f,p_n)=Σ_{k≥0} [Fσ(w_n) S_k(f)−T_k(f)]/w_n^(2k+2).

For k=0,1,2,3 evaluate the functional coefficient row s_k or t_k in the exact AF moment-corrected basis. Its energy-dual norm is sqrt(v* H^(-1) v). To enclose it, interval LDL factorization of H=L diag(p_i) L* gives positive pivots and interval forward substitution y=L^(-1)v, with the upper bound sqrt(Σ_i |y_i|_upper² / (p_i)_lower). The arithmetic grid is AF's directed 10^(-55) rational grid; square-root upper enclosures use integer arithmetic with scale 10^100 to avoid a precision floor in the very small tail sums. All upper bound tests are exact fractions.

For k≥4, ||d||₁≤2||f||₂ (8 modes of physical mass ≥2), hence ||d||₁≤7·10^12 sqrt(q₀[f]), by the rechecked AF margin. The residual geometric series is bounded pointwise by

    2·7·10^12·12^8 sqrt(q₀[f]) (30/w_n^9+360/w_n^10),

since w_n≥3n/4 and 1/(1−(12/w_n)^2)<2. For each integer p≥1 the parity-restricted sum of w_n^(-2p) over n≥N₀ is at most

    (4/3)^(2p) / [2(2p−1)(N₀−2)^(2p−1)],

by the decreasing integral test with spacing two. These bounds give explicit ℓ² control for all k≥1 and for the T₀ term, without subtracting an absolute spectral error.

## 4. The leading tail must exploit the original arithmetic oscillations

Using only |Fσ(w)|≤30w for k=0 is too crude at N₀. Put a_n=2Σκ_h sin(nπ log(h)/4). In each parity, geometric summation and Abel summation give

    |Σ_{n≥N₀, n mod 2=p} cos(nπ s/4)/n²| ≤ 1/[N₀² |sin(π s/4)|]

whenever the denominator is nonzero. For s=0 the bound is 1/[2(N₀−2)]. Expand a_n² with sin u sin v=(cos(u−v)−cos(u+v))/2, using s=log h±log k. Every nonzero denominator is separated from zero by AF's exact trigonometric intervals. The diagonal difference terms s=0 are bounded without division. Multiply the absolute bounds by the exact upper intervals of 2κ_hκ_k and sum all ordered pairs. This yields an upper bound W_A for Σ a_n²/n². Triangle inequality with bσ(w)/w≤8/5 then bounds the leading ℓ² norm of Fσ(w_n)/w_n² by

    (8/5) sqrt(Σ w_n^(-2)) + (4/3) sqrt(W_A).

Its audited approximate upper diagnostic is 0.028006607604, uniformly in parity. This is an oscillatory tail estimate, not an assumption of random arithmetic phases.

Combine this bound with the S₀ dual norm, all other k≤3 dual norms, and the geometric remainder by Minkowski. Divide the resulting raw-cosine norm by sqrt(2) for the physical e_n. Since D_h^(-1)<8I, the relative weighted mixed penalty is at most 4 times the square of the raw-cosine amplitude upper constant. AG verifies strictly, by exact rational comparison, <212/1000 (even) and <206/1000 (odd); both <11/50. Diagnostics: 0.2110494443763355 and 0.20500609581196. Orthogonal parity splits extend the bound to V.

## 5. Genuine restricted infinite-dimensional Schur positivity

For f∈V and y∈Y, weighted Cauchy–Schwarz yields

    |q₀(f,y)| ≤ sqrt((11/50) q₀[f] <D_h y,y>).

Since 11/50<1/4, the actual complete form satisfies, for (f,y)≠(0,0),

    q₀[f+y] > (1/2)(q₀[f]+<D_h y,y>)
              > (1/2)·10^(-25)||f||² + (1/16)||y||².

The orthogonal physical split makes this a true positive infinite-dimensional sector containing the near-null trials, conditional on the previous analytic inputs. On the original smooth pole-null core, take the inverse image of V under projection onto the first 32768 Neumann coefficients: for such a smooth source z, f=P_low z∈V and y=z−f∈Y. Its codimension is at most 32768−14=32754 (at most 16377 per parity), so the conditional negative-index upper bound improves to 32754. This is not a measured index. The argument does not require the intermediate complement of V to be positive or a compactly supported approximation to preserve the high sector.

**Limitation:** the intermediate modes below 32768 and outside V, and their mixed interactions with V and Y, are not controlled. The result does not bound the full low-block Schur operator, does not prove the negative index is zero, and is not a global coercivity theorem. No negative original smooth direction found.

## 6. Reproduction and next gate

Run `python3 canonical_return_t1d_ag_relative_tail.py` with the unchanged `canonical_return_t1d_af_neumann_ritz.py` beside it. `audit_output.txt` is the executed local output. The script recomputes the AF interval matrices, all fourteen positivity pivots, the relative functional norms, the all-pairs prime oscillation bound, the infinite geometric remainder and both strict relative mixed-tail inequalities. Mathematical form identities, interval arithmetic and tail derivations remain subject to independent review; not Lean or remote CI.

**Next AH:** exploit the same divided-difference representation to treat the intermediate sector, particularly by adaptive low-energy whitening and a structured positive high-block inverse rather than expanding the entire low matrix at a fixed absolute tolerance. Seek a rigorous spectral-sign enclosure or a genuine original negative witness. A certified V⊕Y sector is not positivity on the entire original source. Both pole directions beyond the pole-null restriction, support-independent positivity, completed Suzuki–Connes transfer, Phase IV and RH remain open.