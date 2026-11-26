# CR-ALP-055 — Complete E850 finite section; missing-residual targets falsified by global bridges

DRAFT / NONPROMOTING. All continuum identities are inherited classical premises. No Lean certification, complete fixed-collar positivity, arbitrary-support positivity or RH result.

## Global bridge and exact source reconstruction

Replayed CR-ALP-054 saved verifier, all 340 shifted pivots PASS. Retain gamma=72/25, L=gamma I+Q85 P Q85, ell=51/200, u=1101/200, and the frozen CR-ALP-053 rational trial X on G=E600 minus E85. Let D=span(e601,...,e850) intersect H0. The previous J=G plus D excluded one genuine global bridge in each parity. In a parity with pole function b, let c_D=(S850-S600)b, s_D=||c_D||², s_T=||R600 b||², r=s_D/s_T. The bridge h=(c_D-r R600 b)/sqrt(s_D(1-r)) is unit, pole-null, orthogonal to J and Q850. J'=J plus span(h) in each parity is exactly E850 minus E85, with complement Q850 in Q85. It is represented by finite sine coefficients (r c_j below 601 and c_j on 601..850) divided by sqrt(s_D(1-r)), and complete pole coefficient -r/sqrt(s_D(1-r)).

Regenerate inherited signed coefficients at 100-decimal outward integer scale using CR-ALP-054's engine (shift 1024, 24 Bernoulli terms). Compute <h,TV> from genuine sine-pole and pole-pole entries, the latter recovered from pinned anchor diagonal. Compute P h, <b,P h>, <g,P h> for old bridge g, and <h,P h> from the original seven signed shifts. Gamma cross term of L with X is zero by physical orthogonality. New residual row r_h=<h,TV>-<h,P X>. Reuse the pinned CR-ALP-054 R and P banks, check interval overlap of reconstructed K and all old residual entries, append exactly the computed bridge row and projected prime cross column/diagonal. No new numerical Galerkin solve; X remains unchanged. Generator runs both parity sectors through 850; redundant coordinate dimensions 384/385, two discarded pole lines each, physical dimensions 382/383, totaling 765 high directions. E85 plus J' is the full 850-dimensional E850.

## Executed exact tests

Independent standard-library replay computes norm of r_h directly from its interval endpoints. Odd squared norm lies in (0.05859,0.05860), even in (0.05502,0.05503). Since h lies in the missing residual sector of CR-ALP-054, its proposed tolerances 1/250 and 1/60 are false for the frozen trial. This alone does not falsify the global affine or true inverse comparison.

Evaluate affine penalty A=(gamma R_J'^* R_J' - R_J'^* P R_J')/(ell*u) and compressed one-sided Chebyshev inverse polynomials in 24-digit outward integer arithmetic. Odd affine shifted 43 pivots PASS at 1/100. Even affine has exact normalized rational negative witness upper <-0.01049. This falsifies that sufficient comparison on the enlarged projected residual, not full original W or the true inverse. Degrees 4,8,12 all pass at 1/16 in both sectors (255 shifted pivots, plus odd affine 43). Thus K-R_J'^* L_J'^-1 R_J' > 1/16 I; since X lies in J', the compressed residual identity yields positivity of the structured finite block and original W on full E850. Higher-degree compressed moments do not bound the global inverse directly.

Independent replay pins both predecessor and new bank hashes, verifies exact retention of old compressed R/P entries, reconstructs affine penalty from saved G/M, checks integer Bareiss pivots and the rational rejection witness, and independently bounds bridge row norms. It does not independently regenerate analytic entries or polynomial moments. Midpoint eigenvectors are used only to select a rational rejection witness, with exact interval acceptance.

## Conditional next correction

Take Z=q12(L_J') R_J' in the augmented space, with q12(t)=(1-r12(t))/t, r12(t)=T12((gamma-t)/(21/8))/T12(192/175). This exact polynomial correction is defined but not constructed as a frozen new trial in this stage. X'=X+Z. Its corrected low energy K' >= K-R_J'^* L_J'^-1 R_J' > 1/16 I. Complete residual R'=R-LZ decomposes orthogonally into r12(L_J') R_J' and E=Q850(R-LZ). The new exact residual Gram row bounds give ||R_J'||²<1/5 in both parities; eta12²<1/8000. Thus in-band residual squared norm <1/40000. If ||E||<=1/10 in each parity, the full corrected Schur margin exceeds 1/16-(200/51)(1/40000+1/100)>1/50. This is only a conditional target, not an evaluated infinite residual. Neither bridge at 850 is missing from this correction. Need compute actual Q850 leakage of X' before claiming full collar. Nonpositive index <=66 unchanged.

## Replay

From archive root run `python cr_alp_055/check_saved.py`. Generator: `OPENBLAS_NUM_THREADS=1 python cr_alp_055/generate.py`; moment evaluator: `python cr_alp_055/moments.py`. build.py reproduces generator by patching predecessor generate.py (included). Predecessor 054 bank and required 046/049/053 inputs included. The standalone exact verifier checks saved arithmetic, not all analytic premises. No receiving acceptance changes.
