# CR-INT-T1D-BG — Marginal Schur gain and 896-carrier odd certificate

**Status:** local directed-grid spectral audits PASS, conditional negative-index upper bound 343 = 192 inherited even + 151 new odd; complement margin 3/2000. Original fixed interval (-2,2), smooth double Mellin-pole-null source, 24 signed prime powers. Conditional on inherited analytical comparisons and ordinary binary64 exact-integer-dot arithmetic model. No Lean, remote CI, independent mathematical review or new frequency scan.

## Marginal-gain proposal
For regularized effective metric D+rI with r=1/256, and existing block o, candidate block n, Schur residual K=D_nn+rI-D_no(D_oo+rI)^-1D_on, R=(L_n-D_no(D_oo+rI)^-1L_o)Q on targeted orthonormal Q. For a single candidate with positive Schur diagonal k, gain is (8/35) R_j*R_j/k on Q; its trace is (8/35)||R_j||^2/k. Successive scalar Schur eliminations update R and K. Exact formula holds for SPD metric; the implemented ranking uses only floating central checkpoint matrices and is a heuristic, not a verified gain enclosure or a global optimum. New trial validity rests on BF/BD directed-grid audits independently of proposal. From BF's regularized comparator 193 numerical eigenvectors >.985, shortlist 24 low and 40 high distinct carriers by leverage. Greedy marginal gain selects 32; separately compare 12 low +20 high leverage selection from same shortlist (same number, not identical composition). Predicted greedy trace gains first32 144.03964688545443 / all64 178.18960502024572 (~81%); first 10.718672275749343, last .45835566796471283. Not a certificate of actual inertia or of ranking superiority.

## Executed expansion and custody
All seven BF projected checkpoint hashes match; all five existing cross-factor row blocks and error vectors preserved exactly. New 64 carriers give 896 total. Five signed cross rows PASS. Incremental rectangular metric checks 2*832*64+64^2=110592, preserve 692224 old interval entries in each metric, all 49 cells and 24 prime powers. Full max row interval width 73401/10^9. Carrier residual 38016073477/35184372088832 <.001081. Effective cross-factor error 122926937495369841/36893488147419103232 <.003332, effective metric error 224792457471305779/322818021289917153280 <.000697, carrier Schur correction error 205632780162541619/73786976294838206464 <.002787. Full metric checkpoint PASS. Subsets use principal restrictions of same enclosed matrices, retaining full operator error allowances.

## Exact-grid spectral acceptance
All use ridge 1/256 and BD ridge-absorbed error with all base/product allowances retained.

| Trial | Carriers | Absorbed | Odd constraints | beta diagnostic |
|---|---:|---:|---:|---:|
| BF inherited | 832 | 4950 | 193 | .9892853965 |
| BG marginal greedy | 864 | 4971 | 172 | .9916358231 |
| BG stratified leverage | 864 | 4972 | 171 | .9963163218 |
| BG all shortlisted | 896 | 4992 | 151 | .9959311018 |

Exact betas: greedy 3556168871205937972381000839975635649810975645914637/3586164182982095678289775089903319145463535395405824; leverage 3572954086944971406502637592191330348379648715362411/3586164362275157010511007731679588341132326810419200; all 41668335601565087382614355215569359400382184541373621/41838572493657256734389835373224401251278472337162240. All <997/1000. First leverage verification timed out; separate rerun PASS. Leverage outperformed trace-greedy by one dimension at matched size; no superiority claim. Inherit BE even 192 (including constant); combined conditional index <=343, 42 fewer than BF. Physical margin remains 3/2000. Upper bounds, not observed negative eigenvalue counts. No new joint Mellin pole projection.

## Reproduction
Use BF dependencies and actual checkpoint in ../cr_int_t1d_bf; link only read-only input files, never outputs. Run bg_select.py, bg_cross.py F/G0/G1/G2/G3, bg_metric.py cross/new/combine, bg_schur.py, bg_prebuild.py, bg_finish.py, bg_rank.py, bg_baseline.py, then OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=256 python bg_ridge.py, bg_selected_ridge.py, bg_leverage_ridge.py; bg_check.py verifies custody and recorded rational comparisons. Projected archive includes actual full enclosed metric, carriers, new rows, interval blocks, ranking, near-subspace proposal and manifests. Selected and baseline checkpoints are reproducible by restriction. Large upstream matrices needed only for full rebuild, not projected rerun. Original scripts inherited exact integer-dot model; no claim of independent hardware/formal verification.

## Next BH
Mirror 64-carrier batch in even parity. Consider threshold-focused objective (e.g. score near beta rather than total trace) and/or full batch expansion; do not treat greedy trace as proven better. Explicit joint pole-null projection and constant coupling remain open. Fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.
