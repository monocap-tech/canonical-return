# CR-INT-T1D-CD — signed-symbol secant compression

Executed local matched three-even and twelve-odd trials on the actual CA/CB checkpoints, conditional on inherited analytical comparisons, form/core transfer, and directed integer-dot model. No new dimensions. Conditional negative-index upper bound 225 unchanged. No Lean or independent continuum proof review.

For omitted symbol J=mI+D, ||D||<=r, the exact scalar secant inequality J² <= 2mJ+(r²-m²)I gives, for w=v+h and s=m+1/2,
<u,(J+J²)u>-2 Re<Ju,w> <= (m²+m+r²)||u||²-2m Re<u,w>+2 Re<Du,su-w>.
The last term is bounded by r(theta||w-su||²+theta^-1||u||²). The full shifted Gram is Gv+Gh+s²Gu-s(Gvu+Gvu^T) because physical interior/exterior supports are disjoint. Subtract the retained Gram of Z+Y-sU with directed rational error. The coefficient m²+m+r² is checked nonnegative before replacing the omitted u Gram by its upper envelope. Exact fractional centers are divided without integer truncation. Keep the CC directional action and mixed-energy allowances. Heuristic float grid only chooses theta/gamma; fixed replay uses rational LDL, no eigensolver.

| sector | CC directional | CD exact normalized | strict rational margin |
|---|---:|---:|---:|
| even 3 | .3102199999994317 | .32874999999939775 | 16437/50000 |
| odd 12 | .04670999999970831 | .049269999999692324 | 2463/50000 |

Theta/gamma even 1/15,1/65; odd 1/41,1/175. Fixed replays 2/2 PASS, exact saved fractions and all shifted grid matrix entries reproduced. The shifted Gram row-upper diagnostics are 8.472225346140151e-6 even, 5.634798602268931e-5 odd; compare CC's different unshifted Gram, not an identical operator. No new source/action generation. Reuse the six pinned archives in custody/cc/cc-inputs.sha256 and the original extraction layout described in CC. Copy cd_secant.py and cd_results/ into common parent; run OPENBLAS_NUM_THREADS=4 python cd_secant.py even --fixed and odd --fixed. Exact independent terminal checks verify the positive principal minors only, not source generation.

## Four-column even frontier

Actual CA saved proposal has Q(121,8), X(48,121), v0(121,8), g0(1,8), t(1890,8), T(2048,8). Numerical central values descending .94036453,.59644899,.44435553,.00282928,-.00270105,-.00704257,-.03417233,-.07315323. Fourth direction selected but not enclosed. Need three BK inputs, not two:
- even/exceptional.npz SHA256 4c5c5b984027b93748450b03b04fdeeab27762a71844554d8d963e992008e7f9
- even/physical_dictionary.npz SHA256 bf606cf0cd9f9c64ca8d8148a468f77706c7e9c866c128775627cefba019749b
- even/pole.npz SHA256 b890700e1ed8103df55ab857562843c64b4f05fbd565ecb113497b75abb4acb2

These hashes are from CA verified_input_sha256.txt. At archive custody commit 42ff1bd357a5e28e18fc3d98d3950dac7811bdaf, the BK four chunks are present and total 2,992,441 bytes. Read back last binary chunk in base64 and parsed actual ZIP central directory: all three required member names confirmed, uncompressed sizes 1,335,567; 33,555; 304,250 bytes respectively. This is not reconstructed member hash verification. Library materialization still denied; local GitHub DNS unavailable. Full fourth-column generation not executed. On a separate scratch copy, verify all three member hashes, change all four [:,:3] selections in CA bx_prepare.py to [:,:4] (including T selection; count must be five), run bu_circle.py 0,1,2, bx_combine.py, bv_joint_action.py even 2048, then use CD comparison on generated data, preserving original manifests and separately hashing regenerated files. Do not assert positivity on a failed lower comparison. No changes to the active archive upload branch. No global positivity or RH conclusion.
