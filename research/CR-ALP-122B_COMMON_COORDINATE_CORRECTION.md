# CR-ALP-122B — Common-coordinate and correction transport audit

Research-only, nonpromoting. CR-ALP-120 establishes conditional rank 85. This note separates executable source identities from inherited continuum premises. The remaining historical archive generation must be reopened; matching a saved matrix hash is not proof of its analytic identity.

## Source trace (reopened actual code)

- 074/075 source bank (6a644138750ba8d4e186a8a5878c4bbee97cf022, 61110f405a5d7339e076dbdfea86e2d984169ed7): pin the three 057 source/group/trial banks, use the same zero-based coordinate column for all source coefficients, pole coefficients, trial D/Db and grouped boundary rows; output coordinate-indexed hashes and omitted group norm upper bound.
- 082 inertia receiver (d478188f8209227e2d7a8f9eb638594d50d88e3a): reconstructs inverse bounds from coordinate-labelled logs; requires all 85 source hashes; requires complete coordinate orders J as permutations.
- 085 obstruction block (e450dde6a33520f13c0edaef1228b5168daaf0b2) and 086 terminal moment (2a9d96ec90d6b71bd9dc3c0b786af0636a0c6e32): same permutation J for both K axes and physical input file order; combined inputs use signed interval linear combinations of all coefficients; source errors scale with absolute coefficients. Subsequent rational directions are not a new physical basis.
- 075 shifted energy (d8958b1f457f711b1ac14562fb6335752cc6284f): intersects transposed K intervals, forms midpoint grid and subtracts row-radius sums. Hence the lower matrix is a Loewner lower bound of the original saved K in order J, conditional on entry enclosures.
- 091 coupled capped (0bfa6d7761d7153bb3da976dfbe8468e816107c1): reopens the same 057 trial bank; obtains tau from the complete physical high-pole norm and finite sine precursor. The exact correction Gram includes the infinite pole tail via Parseval; no zero-tail substitution.
- 092 receive block (dab163b85105f2c4fb0fc5e24433cab458533ddc): same K lower, correction Gram and obstruction coordinates, and a single polarized quadratic majorant.
- 053 historical note identifies V as the frozen 049 energy coordinates, and defines K=V^*H V-B_V^*X-X^*B_V+X^*L X. 055 and 056 retain V while enlarging and correcting X. This is historical provenance, not a direct verification of the archived generator.

## Exact correction identity

For each parity let S be the first 85 sine projection in that parity, h=(I-S)b and st=||h||^2. On the constrained high space, Q=(I-S)-|h><h|/st. For a finite precursor t=sum_j D_j e_j (j<=850), X=Qt, hence X_j=0 for j<=85 and X_j=t_j-tau c_j for j>85, tau=<h,t>/st. All coefficients are real for the archived coordinate columns. The saved Db overlap with -tau is a necessary consistency check, not sufficient proof of exact generator equality.

Let A_j=min(Dcap,D_j^Gamma)-gamma; A_j=Acap for j>=2049. For a vector of coordinate columns, write row t_j, row tau, and c_j=<e_j,b>. Then exactly
X^*AX=sum_{86<=j<=2048} A_j (t_j-c_j tau)^*(t_j-c_j tau)
+Acap (st-sum_{86<=j<=2048}|c_j|^2) tau^*tau.
The sums use the matching parity only. This is the 091 weighted-row Gram with its final Parseval row; it includes all j>=2049, not just measured modes. Given the pinned coefficient and pole enclosures, outward entry/row reserves prove GL<=X^*AX<=GU.

## Common comparator transport

For L>0, A>=0, C=L+A on Q, define r_L=R-LX, r_C=r_L-AX (A compressed to Q). For arbitrary z and y, 2Re<r_C z,y>-<y,Cy> <= <r_L z,L^-1 r_L z>+<Xz,AXz> by completing the L and A squares separately, without inverting A. Supremizing y gives B_C<=B_L+X^*AX. This is the exact Loewner justification of inherited global and coordinate bounds. It does not require equality of the two response Grams.

For any rational source coordinate vector z in order J, the physical input constructor uses exactly the coefficients of z on the permuted source columns. Thus each valid bound produced on that input concerns the same compression U^*B_C U, provided the original 049/056/057 generator attachment and physical validation hold. A real symmetric certificate extends to complex sources by real/imaginary splitting, using 122A parity and reality.

## Remaining explicit blocker

Reopen the actual authenticated CR-ALP-058 archive (reassembler pins 62,872,042 bytes, 90 chunks, SHA256 f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2). The 053/055/056/057 historical notes alone cannot establish the 049 coordinate map's exact invertibility, the precise 056 K/057 trial generator linkage, or all numerical enclosure assumptions. Stage a source-only extraction and manifest census; do not rerun physical integration or infer generator identity from filenames. Tail Gamma analytic replay is still open. No fixed-collar promotion.
