# CR-INT-T1D-BF — Adaptive odd carrier expansion and incremental metric

**Status:** local directed-grid spectral audits PASS, conditional negative-index upper bound 385 = 192 inherited even + 193 new odd, physical complement margin 3/2000. Original fixed interval (-2,2), smooth double Mellin-pole-null source, 24 signed prime powers; conditional on inherited analytical comparisons and ordinary binary64 exact-integer-dot arithmetic model. No Lean, remote CI, independent mathematical review or new frequency scan.

## Checkpoint custody and proposal
All 17 AZ/BC derived checkpoint digests match. From the regularized 640-carrier BD odd comparator at ridge 1/128, 354 numerical eigenvectors with eigenvalues > .985 give a leverage proposal of 64 additional low and 128 high carriers. Distinct from the previous 640; 832 total. Carrier selection is heuristic, not a certificate of negative directions. All five inherited arithmetic cross row blocks and error vectors match exactly; added row/column errors accumulated as in BE. First combined cross/solve run timed out; separate solve completed. True carrier residual 36964539843/35184372088832 < .001051. Effective cross-factor error 119728579063401783/36893488147419103232 < .003246; Schur correction error 195355839627222787/73786976294838206464 < .002648.

## Incremental exact quadratic metric
For c=8/35 and J0=(399/50 I-A)/4, K=P*J0P/21+c(J0P)*(J0P). Preserve the entire old 640x640 interval block (409600 entries), rather than recompute it. The full square integration timed out; bf_metric.py adapts BB's cell antiderivatives to rectangular pairs. Evaluate old-new and new-old separately (122880 entries each) and new-new (36864 entries). All 282624 new AX overlap checks pass across the 49 cells and all 24 prime powers. Intersect transpose enclosures and assemble symmetric intervals, verifying the old block unchanged. Full max row interval width 34019/500000000. Effective D=K-cQ_Y enclosure error 213091649378017027/322818021289917153280 < .000661. No pilot-relative metric percentages claimed.

## Executed spectral audits
Use BD's ridge-absorbed error and exact-grid Gram/residual test with m=5143, n=12811, 832 carriers. All base and directed product error terms retained.

| Ridge | Absorbed odd dimensions | Remaining constraints | beta diagnostic |
|---|---:|---:|---:|
| 1/128 | 4939 | 204 | .9952320444040321 |
| 1/256 | 4950 | 193 | .9892853965316815 |

Exact accepted beta at ridge 1/128 = 417516243702352262679897660362334044219875171367909741/419516479649095930116715444412451193462695922132582400; ridge 1/256 = 1609913160097440788666529260535678353684391647828457/1627349565394988592959915305855949620574216058306560. Both <997/1000. The initial 4951-dimensional proposal at 1/256 failed (.9977890670615082); nested 4950-dimensional test passed. Combining inherited BE even 192 and new odd 193 gives conditional n_-(q0;D^0)<=385, 160 fewer than BE. q0[f]>=3/2000||f||^2 on certified complement. These are upper bounds, not observed negative eigenvalue counts. No new even execution or joint pole projection.

## Reproduction
Use BC's odd helper scripts and arrays and actual checkpoint, plus BF scripts. In separate writable BF directory link/copy read-only upstream derived matrices (verify manifests), do not link output files. Run bf_select.py, bf_cross.py F/G0/G1/G2/G3, bf_schur.py, bf_prebuild.py, bf_metric.py cross/new/combine, bf_finish.py, bf_ridge.py with OPENBLAS_NUM_THREADS=4 and RIDGE_DENOM=128 or 256, bf_check.py. The parent path in bf_metric.py and bf_check.py is ../cr_int_t1d_bc; bf_cross.py and bf_select.py support BF_PARENT. Projected archive contains actual enclosed effective metric, metric intervals, selected carriers and incremental cross rows. Extract that archive plus BF scripts and inherited helper scripts/source arrays into a separate working directory to rerun bf_ridge.py. Original upstream matrices remain locally available but are not part of projected archive; do not assume future local persistence. No global optimum for carrier selection or ridge claimed.

## Next BG
Both parities now have 832 carriers. Test exact marginal Schur gain rather than raw leverage, especially on remaining 192/193 constraints, and preserve rectangular metric updates. Incorporate original Mellin pole-null projection jointly (including even constant), not simply subtract a constraint count. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.
