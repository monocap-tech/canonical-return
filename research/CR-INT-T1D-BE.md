# CR-INT-T1D-BE — Adaptive arithmetic-carrier expansion

**Status:** even expanded spectral audit PASS, odd BD inherited. Conditional negative index <=545 (192 even + 353 odd), physical complement margin 3/2000. Original fixed window (-2,2), smooth double Mellin-pole-null source, 24 signed prime powers. Inherited analytical comparisons and ordinary IEEE binary64 exact-integer-dot arithmetic model. No full frequency rescan, Lean, remote CI, or independent mathematical review.

## Carrier expansion and custody

All 16 BD derived checkpoint hashes pass. From the approximate regularized BD comparator at ridge 1/128, select the 356 numerical eigenvectors with eigenvalues >.985. Coordinate leverage sum of squared components proposes 64 additional low and 128 high deficit carriers, disjoint from the 640 old carriers. Keep old row order and append new carriers: 832 total. Selection is heuristic, not evidence of negativity or global optimality. The variational family on the enlarged carrier range contains the old one by zero-padding the trial coefficients; the exact optimized upper bound cannot increase. Computed certificates need separate error checks.

Only new arithmetic rows are evaluated: F and G0–G3, with existing exact directed phase tables. Combined row errors concatenate and maximum column errors add. All five old row blocks, row error vectors and input IDs agree exactly with BD. Save the added rows and bounds in be_cross_delta.npz. The first metric integration timed out; restart completed. Full 832-carrier piecewise metric across 49 cells and all 24 prime powers passes 692224 independent AX overlap checks, maximum interval row width 17413/250000000. Gershgorin relative-metric pilot percentages are not established. The full 832-RHS solve passes with residual 421550049/549755813888 < .000767 and trial norm <=6926089/2097152. Split prebuild and finish checkpoints establish eta_L=1164887664407803/576460752303423488 < .002021, eta_Q=2177251784279529/1152921504606846976 < .001889, eta_D=2462416137897449/5044031582654955520 < .000489. Inherit BD's ridge-absorbed error inequality, full base enclosure, directed product and spectral residual checks.

## Executed spectral acceptance

| Ridge | Absorbed even dimensions | Even constraints incl constant | beta diagnostic |
|---|---:|---:|---:|
| 1/128 | 4873 | 208 | .9961026394580837 |
| 1/256 | 4889 | 192 | .995683939304389 |

Both runs exit zero. Exact strongest beta=155952017581908153577771242178418805244876153917963/156628033681913499921754309444722091089173501968384 <997/1000. Thus 5080-4889+1=192 even constraints, and inherit BD's 353 odd constraints. Conditional index <=545, improvement 162 from 707. Coercive margin >=3/2000 on the new complement. Not an observed count of negative eigenvalues.

## Pole-projected constant (scalar certificate only)

The even constant u0=1/sqrt(8) on I, g=cosh(x/2), can be replaced on the pole-null source by u0-proj_g u0. Its squared norm is exactly 1/2-2 sinh(1)^2/(2+sinh(2)). Inherit b0=2041616723/50000000 from AW's constant diagonal lower bound. Exact rational Taylor bounds in be_pole.py establish b0||u0-proj_g u0||^2<2/5 (diagnostic .37179254727122724). The original constant alone would have deficit norm squared b0/2. Mixed coupling with other deficits remains unverified; the constant is still counted among the 192 even constraints. General constrained inverse identity for the positive physical D: T*D^-1 T - (T*D^-1 g)(g*D^-1 g)^-1(g*D^-1 T) gives the exact Gram on the pole-null subspace. Can equivalently replace T by T-g a with any row a in a variational upper bound on the constrained source. Neither mixed certificate nor additional count claimed.

## Next BF

Mirror the targeted carrier expansion in the odd sector; consider a more precise marginal-gain selection. For positive regularized D partitioned into old/new carrier blocks, the exact extra optimized gain before multiplying c=8/35 is (L_new-D_no D_oo^-1 L_old)* (D_nn-D_no D_oo^-1 D_on)^-1 (L_new-D_no D_oo^-1 L_old). This follows by block completion of squares and is a proposed selection heuristic, not an executed new certificate. Test pole-projected constant jointly, and actual constrained inverse for both parities. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.

## Reproduction

Projected archive stores actual new effective matrix, metric intervals, carrier IDs, incremental cross rows, manifests. Use BD reproducibility bundle's helper scripts and arrays with the BE scripts. To run projected gate only, extract projected archive and run OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=256 python be_ridge.py. To regenerate upstream, reconstruct BD's full derived checkpoints as documented there (verify BD manifest), copy/symlink read-only input files into separate BE directory, set BE_PARENT to absolute BD directory; do not symlink output paths. Run be_select.py, be_cross.py F and G0–G3 separately, bb_metric.py az_carr.npy, be_schur.py, be_prebuild.py, be_finish.py, be_ridge.py at both ridges, be_pole.py and be_check.py. The added environment-variable path support was syntax checked after the actual run. Old full-size derived matrices and new Schur checkpoint remain in local working directories, not included in projected archive; cross-session persistence not guaranteed. Archived appended rows allow reconstruction without recomputing phases. Full upstream scripts are inherited from BD, no claim of revalidation of original scan.