# CROSS-PW-004 — moment-compatible fixed detector and terminal supplier handoff

2026-11-26 PST. Mathematical detector gate: classical PASS. Native attachment: OPEN. Comparative lane: CLOSE / HANDOFF (not a proof of RH, not a completed window supplier). CROSS-SD direct substitution remains closed. No new Lean declarations, external code imports, or explicitly dispatched workflows; no mutation of receiving branches.

## Source and custody dependencies
Stage 003 canonical-return@b09ee4084a900031b2edeece2372d77f30822ddd research/CROSS-PW-SD-003_LOCALIZED_MIXED_RECEIVER_AND_HANDOFF.md blob 36332bcefd99cfa0b7e0a748fb4a2c5a627d50fc, directly reopened. Stage 001/002 pins retained. RETURN106D original separated form, at 6c0f6eeb133a0c5ee3d7420c2861c693deea4006 blob 045c9afa3914832dd3510f734f3bd5a322a62a5a. General Landau source leanprover-community/riemann-hypothesis@40e56bab026552e089d7afcfeff08d2defb9c005, src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean blob 7915058f10b33a52652290cd887cbccdc00e676a (literal theorem reopened in interrupted stage 004); PoleData and WindowHypU source pin in 003. RECON weil-lab@3973964c66f368d1aa3fbe58dda833245e879dfe WeilDefect/MasterRecon/V10B52b.lean blob 4a2e6b845d9cb1420cc13d4b7e37e244a43392f4, literal original_moment_kernel at t=1/4 and j<15 (source inspection only, no current native audit verified). V10B52a.lean blob afe76f8a911789cb4f5e1e4f92385884dc95639b at same commit. Earlier V08I logCarrier and V07A Mellin translation pins in 002. The construction below is independent of the earlier external universal-sinc detector and does not require its code or identity.

## 1. Fixed smooth compact detector with no off-imaginary-axis transform zeros

For j>=1 set r_j=2^(-j-3), let mu_j be normalized Lebesgue measure on [-r_j,r_j], and define mu as the law of the sum of independent X_j with laws mu_j on the countable product probability space. Since |X_j|<=r_j and sum r_j=1/8, the series converges uniformly on that product of intervals; mu is a real even probability measure supported in [-1/8,1/8]. For all complex w, dominated convergence gives
 M_mu(w)=integral exp(wx)dmu(x)=product_j sinh(r_j w)/(r_j w)
 with factor value 1 at zero. Each factor is entire and its Taylor expansion differs from 1 by O_K(r_j²) uniformly for w in a compact K. Since sum r_j²=1/192, the product converges locally uniformly to an entire function. If Re w!=0, every factor is nonzero (sinh z vanishes only at integer multiples of pi i), and the convergent product is nonzero: in the tail choose factors within 1/2 of 1, use absolutely summable principal logarithms there, then multiply by the finite nonzero prefix.

At w=i y each factor is sinc(r_j y), modulus <=1. For every N>=1 and y!=0, |M_mu(i y)|<= product_{j<=N} 1/r_j * |y|^(-N); combining with the bound 1 gives rapid decay with a constant depending on N. With the Fourier convention exp(-2pi i xi x), the characteristic function is M_mu(-2pi i xi). It and every polynomial multiple are integrable. Fourier inversion for finite measures gives a continuous density h, and differentiation under the inverse Fourier integral of all orders gives h in C-infinity. Since mu is supported in the specified compact interval, h vanishes outside that interval by continuity and equality of measures; h is real even nonnegative and integral h=1. Its bilateral Laplace transform equals M_mu. This uses standard classical Fourier inversion and countable product probability, not RH or estimates on actual zero heights.

## 2. Simultaneous moment and pole annihilation

Let D=d/dx, f=D^16(D²-1/4)h. Then f is real, even, C-infinity, supported in [-1/8,1/8] (strictly inside RECON's [-1/4,1/4] and RETURN's [-5/4,5/4]). Compact-support integration by parts gives M(D^k h,w)=(-w)^k M(h,w). Hence
 M(f,w)=w^16(w²-1/4) product_j sinh(r_j w)/(r_j w).
 Nonzero f, e.g. M(f,1/4)!=0. All monomial moments integral x^k f(x)dx vanish for 0<=k<=15, by differentiating the transform at zero, or 16 integrations by parts against x^k. Both elementary moments M(f,+1/2)=M(f,-1/2)=0. This is one extra monomial moment beyond the 15 constraints j<15 in V10B52b; every polynomial of degree <15 pairs to zero, including the specified Legendre rows. Because f is smooth compact it has finite logarithmic energy with weight log(2+|2pi xi|). An explicit spectral witness for V08I is sqrt(log weight) Fourier(f), square integrable; inverse-weight then inverse-Fourier yields f almost everywhere, so support constraint holds. These are classical carrier identities, not a compiled construction of that Lean subtype.

For every nontrivial off-critical-line zero rho, w=rho-1/2 has Re w!=0, w!=0, w!=+-1/2 because 0<Re rho<1. Thus M(f,w)!=0. Evenness gives M(f,-w)=M(f,w), so the coefficient m_rho M(f,w)^2 in the separated mixed explicit formula is nonzero (positive multiplicity m_rho). A single fixed detector sees all off-line zeros without choosing its scale after the zero or requiring simultaneous windows for a family. No claim about nonvanishing at on-line zeros, which need not be detected.

For any real translation a, M(tau_a f,w)=exp(aw)M(f,w), so the two pole constraints and zero of order >=16 at w=0 remain unchanged. Equivalently expand (x+a)^k by the binomial theorem. The logarithmic norm is translation invariant by the unimodular real-frequency Fourier phase. Translating a test does not preserve a fixed centered support aperture; use the enlarged actual carrier.

## 3. Exact final supplier in a moment-constrained two-lobe sector

Use Stage 003 direct Landau receiver with this f. Let U>2*(1/8), take any fixed delta>0, e.g. 1/2. Define
 v_u=(tau_{-u/2}f-tau_{u/2}f)/sqrt(2), u>1/4.
 It is a real odd smooth test, supported in [-(u/2+1/8),u/2+1/8], all moments k<=15 and both elementary pole moments zero. Disjoint lobes give ||v_u||_2²=||f||_2². Translation invariance and triangle inequality give ||v_u||_log²<=2||f||_log², uniformly in u. The original real symmetric Weil form is translation invariant, so exactly
 W(v_u,v_u)=W(f,f)-W(f,tau_u f).
 Put [x]_+=max(x,0). The positive-part map is 1-Lipschitz, hence
 |[W(f,tau_u f)]_+ - [-W(v_u,v_u)]_+|<=|W(f,f)|.
 Thus the two weighted integrals over [U,(1+delta)U] differ by at most |W(f,f)| exp(-epsilon U)/epsilon, for epsilon>0. Allowing eta/2 in one receiver and choosing U sufficiently large absorbs the fixed additive allowance into exp(eta U). Therefore the Stage 003 RH-equivalent receiver is also exactly equivalent to:
 for every 0<epsilon<1/2, eta>0, and arbitrarily late U,
 integral_U^(3U/2) exp(-epsilon u) [-W(v_u,v_u)]_+ du <= exp(eta U).
 Here a cofinal integrated lower-energy allowance on one fixed odd two-lobe family suffices, rather than full arbitrary-support positivity or two-sided control of every cross block. This uses the full original Weil form and Stage 003 meromorphic receiver. Conversely RH implies the condition; that implication is NOT a permissible prospective supplier. The minimum separation in the original RETURN106D note uses a=5/4; simply require U>5/2 when applying that source verbatim.

The actual window budget remains OPEN. Fixed-support moment-kernel equivalences, local Schur certificates, bounded sampling with aperture-dependent constants and finite shift censuses do not establish it. The construction is not a claim that a negative witness exists. No numerical evaluation of actual detector energy or large-scale prime discrepancy is asserted. No arbitrary-support positivity result.

## 4. Terminal acceptance and handoff

PW4-D: classical detector/constraint/visibility construction PASS. Native admission, translation, all-zero visibility, full signed explicit formula and general Landau attachment still OPEN; do not count as accepted Lean declarations. PW4-S: actual cofinal weighted negative-part budget on literal v_u OPEN. Assign conceptually to RETURN (arithmetic interaction supplier) and RECON (moment-constrained spectral/energy supplier), reuse rather than start an independent eighth lane. No actual branch edit, issue assignment, CI dispatch or automatic worker initiation performed. Stop the comparative investigation; reopen only for concrete independent supplier or native transport work authorized in a receiving lineage. SD direct substitution remains CLOSE, moment invariants archived. No RH proof claimed.

## 5. Reproduction
verification/cross_pw_004.mjs commit 0448dc0c06c16dba9ddb3b78c6d3ec7bcdeddf2b, blob 80078992f4287c7a9faf6cfc4a22b2f2c58049c8. Executed the identical replay function with BigInt arithmetic: 128/128 radius partial sums, 16/16 moments of finite polynomial model 4D^18(1-x²)^20-D^16(1-x²)^20 on [-1,1], 2/2 pole symbol checks, nonzero even model PASS. The finite polynomial is not h or f above; these checks corroborate finite algebra only. The infinite convolution, smoothness and zero-visibility proof is the classical derivation in this record. No native Lean or numerical prime calculation. Exact remote script readback passed. Verify this note and head before closing custody.
