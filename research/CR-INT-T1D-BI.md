# CR-INT-T1D-BI — Threshold-weighted bilateral 960-carrier expansion

**Status:** executed local directed-grid audits PASS. Conditional negative-index upper bound 246 = 127 even + 119 odd, physical complement margin 3/2000. Original fixed interval (-2,2), smooth double Mellin-pole-null source, all 24 signed prime powers. Conditional on inherited analytical comparisons and ordinary binary64 exact-integer-dot arithmetic model. No new frequency scan, Lean, remote CI or independent mathematical review.

## Odd ridge refinement
All ten BG projected digests match. Rerun BG's existing ridge audit at 1/512. Initial 4998 proposal rejected (beta .9975962853518173); nested 4997 passes, 146 odd constraints instead of 151. Exact beta 14676661722158924764736602507825939140390825819523209/14742999136319513163904118087507160933587501264666624 < 997/1000. Combined with BH even 156 gives intermediate bound 302. No new carriers in this test.

## Bilateral threshold-weighted shortlist and expansion
Use current 896-carrier parent central regularized comparator at ridge 1/512. Compute numerical eigenvectors Q with eigenvalues v > .985 (156 even, 147 odd). For each unused coordinate j rank by sum_l |Q_jl|^2/(1+((v_l-.99)/.05)^2). Add 24 low and 40 high carriers per parity, disjoint from old 896. This is threshold-weighted leverage, not a certified secular crossing predictor or the full greedy marginal Schur algorithm. Entire 64-carrier batches are tested. Matched-size secular comparison not executed. All 20 parent checkpoint digests match, all five old cross-factor row blocks and error vectors preserved per parity, and all 802816 old metric interval entries preserved per parity for D/Q/K. Incremental metric directed enclosure across all 49 cells and all 24 prime powers passes 118784 new entry-overlap checks in each parity (2*896*64+64^2). Row interval widths 80953/10^9 even, 78849/10^9 odd.

| Quantity | Even | Odd |
|---|---:|---:|
| Deficit dimension | 5080 | 5143 |
| Surplus dimension | 12619 | 12811 |
| Carriers | 960 | 960 |
| Carrier residual | 222923187/274877906944 | 9696263259/8796093022208 |
| Cross-factor error | 610066516994721/288230376151711744 | 31235892795860335/9223372036854775808 |
| Metric error | 675318983163281/1261007895663738880 | 14530949983303923/20176126330619822080 |
| Carrier Schur error | 592734815124881/288230376151711744 | 13245079724554483/4611686018427387904 |

Even combined solve/prebuild/finish/ridge execution hit time limit after writing completed effective metric; separate spectral execution passed. Odd stages completed separately.

## Executed spectral acceptance
Both use BD ridge-absorbed error, ridge 1/512, exact-grid Gram/residual checks; all base and directed product error terms retained. Even absorbs 4954, leaving 5080-4954+1=127 including constant. Exact beta 104723612018048935527928997373268864277694856579/105726999423929531292502718746319934745550520320 (.9905096388685226). Odd initial 5025 proposal rejected (.9976705509016175); nested 5024 passes, leaving 5143-5024=119. Exact beta 113175128818808878723413990314379160309191198259719/114423289108811545293748368466542411603078912409600 (.9890917286181512). Both exact betas <997/1000. Conditional index <=246, 61 fewer than BH's 307, physical margin >=3/2000. Upper bounds, not observed negative eigenvalue counts. Constant remains excluded; no joint pole projection.

## Reproduction and custody
Use upstream BH/BG actual checkpoints and inherited helper scripts and source arrays. setup.py adapts source scripts to separate writable even/odd directories (run in /mnt/data with parents present; never overwrite upstream matrices). For each parity, run bi_select.py; bi_cross.py F/G0/G1/G2/G3; bi_metric.py new/cross/combine; bi_schur.py; bi_prebuild.py; bi_finish.py; OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=512 python bi_ridge.py. Set OPENBLAS_NUM_THREADS=4 for all scripts. bi_check.py verifies custody, preservation, exact recorded comparisons, and writes new manifests and cross deltas. Projected archives include actual enclosed effective matrices and metric intervals, incremental rows, carrier IDs and near-space proposals. Extract corresponding archive with scripts and inherited helpers into separate directory to rerun spectral gate. To reproduce old odd ridge, run odd_old_ridge.py in original BG directory with ridge 512. Prepared bi_threshold.py scripts are unexecuted; no secular-selection result claimed. Full large upstream inverse matrices needed only for rebuild, not projected rerun. No guarantee of future local persistence.

## Next BJ
Compare matched-composition secular/threshold-focused and threshold-weighted leverage selections on the actual 960-carrier checkpoints; consider another 64-carrier expansion or joint pole-null inverse correction. At 246 constraints, direct signed physical exceptional-sector Schur control becomes increasingly relevant. Do not equate declining upper bounds with measured negative eigenvalue counts. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.
