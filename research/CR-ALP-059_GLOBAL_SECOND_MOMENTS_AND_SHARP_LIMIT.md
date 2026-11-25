# CR-ALP-059 — native second moments and sharp scalar limitation

DRAFT / NONPROMOTING. Native verification run 45546466147 at commit 4b60bff3e0bde60d186b216741f5a2f68e15c718 passed both sector integration jobs 138977181118 and 138977181479 and endpoint job 138978239959. Archive reassembly and predecessor replay passed. The native jobs reported an exact 898-panel half-domain mesh. The extracted five-line stdout records and exact endpoint JSON are preserved in src/audits/cr_alp_059/evidence; each was fetched back and verified. These are conditional physical-space enclosures under the retained 058 analytic and forward-error premises, not Lean results or original Weil negative directions.

Actual N=||P_Q r||²: odd [0.29588696822533395,0.2994787991614581], even [0.2833380380256938,0.2861475803777602] (decimal approximations of exact rational endpoints). Both fixed quadratic sufficient tests remain undecided. Native inverse bounds: odd [0.05751140189536849,0.18966557014931398], even [0.05179752412072366,0.19033131508443135]. Available corrected energies approximately 0.121547573369 and 0.084343420038.

## Optimal scalar limitation

Write m1=gamma G+M, m2=gamma² G+2 gamma M+N, with spectral interval [ell,u] as in 058. The optimal upper and lower scalar inverse bounds over positive measures with these exact three moments are

U=G/ell-(m1-ell G)^2/[ell(m2-ell m1)],
V=G/u+(u G-m1)^2/[u(u m1-m2)].

They follow from the quadratic reciprocal majorant/minorant, optimized in its interpolation parameter. Sharpness follows from the two-point measures on {ell,(m2-ell m1)/(m1-ell G)} and {(u m1-m2)/(u G-m1),u}: for an admissible exact triple these nodes lie in [ell,u] and the weights obtained by matching G,m1 are nonnegative by m2 G>=m1². The nondegenerate denominators are verified by the exact interval audit.

Using exact rational input endpoints, sharp lower V is enclosed in [0.06723958,0.06747339] odd and [0.06157090,0.06183551] even; sharp upper U in [0.17958825,0.18186700] odd and [0.16834265,0.17110913] even (decimal display). In both cases V_upper<K_lower<=K_upper<U_lower. Thus for every admissible exact triple in the measured box, there are abstract spectral measures with that same triple giving opposite comparison signs. Refining the existing scalar moments alone cannot close these witnesses. This does not establish the actual P_Q spectral measures equal any extremizer, or that the original Weil form is negative.

## Conditional cubic targets

Let T=<P_Q r,P_Q P_Q r>=<r,P_Q^3 r>, H(t)=(t-ell)(u-t), gamma=(ell+u)/2, s=(u-ell)/2. For a>0,

p3(t)=(ell+u-t)/(ell u)+H(t)(t/a²-2/a)/(ell u),

p3(t)-1/t=H(t)(t-a)^2/(ell u a² t)>=0 on [ell,u].

Hence writing h=s²G-N and f=gamma G-M-ell u K, the cubic upper bound lies strictly below K when T>gamma h+s²M+f a²-2ha. Rational a=2.517 (odd) and 2.134 (even) give sufficient acceptance thresholds T>=0.073 and T>=0.269 respectively, after checking every interval corner. Separately x=alpha r+beta P_Q r gives the lower inverse bound 2<r,x>-<x,Lx>; rational (alpha,beta)=(0.671,-0.589) and (0.542,-0.310) give sufficient comparison rejection T<=-0.568 and T<=-0.383 respectively. The exact thresholds are certified by sharp_scalar_postcheck.py and its result JSON. T is NOT measured.

Replay from repository root: python3 src/audits/cr_alp_059/sharp_scalar_postcheck.py src/audits/cr_alp_058_recovery/verification_result.json src/audits/cr_alp_059/evidence/second_result.json /tmp/sharp_scalar_result.json ; compare with retained evidence/sharp_scalar_result.json. Local Fraction execution used exact input entries copied from fetched 058 JSON and native 059 endpoint log; a fresh checkout replay of this new postcheck has not yet been submitted.

A parallel CR-ALP-060 matrix-moment candidate (commit 95f9af940362b24ab41b7c9c7528eb0681b70494, research/cr-alp-060-matrix-moment-bound) derives a full matrix upper/lower comparison and includes exact abstract regression tests. It can exploit cross-column information unavailable to three scalar moments. It has not evaluated the full bank. Next scalar step: reconstruct v=Q85 P r in physical space and integrate <v,Pv> with rigorously transported errors, before deciding whether a full matrix moment bank is useful.

Custody: native outputs and endpoint JSON, new postcheck script/result and this note committed to the original draft branch; immutable 058 archive retained. No promotion, no mathematical claim beyond this comparison. Historical standalone archive reconciliation remains open.
