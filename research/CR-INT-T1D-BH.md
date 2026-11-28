# CR-INT-T1D-BH — Even expansion and threshold-focused carrier selection

**Status:** executed local directed-grid audits PASS. Conditional negative-index upper bound 307 = 156 new even + 151 inherited odd, complement margin 3/2000. Original fixed interval (-2,2), smooth double Mellin-pole-null source, all 24 signed prime powers. Inherited analytical comparisons and ordinary binary64 exact-integer-dot arithmetic model. No Lean, remote CI, independent mathematical review or new frequency scan.

## Full even expansion
All four BE projected checkpoint hashes match. The previous regularized comparator at ridge 1/256 proposes 193 numerical eigenvectors above .985. Add 24 low and 40 high carriers by leverage, distinct from existing 832. All five inherited cross-factor blocks and error vectors and all 692224 old metric entries are preserved exactly. Rectangular directed integration of all 24 prime powers over 49 cells passes 110592 additional AX entry-overlap checks; row width 75377/10^9. Carrier solve residual 54345177/68719476736 < .000791; effective cross-factor error 149346287164451/72057594037927936 < .002073; metric error 40475736092569/78812993478983680 < .000514; carrier Schur correction error 35664198315929/18014398509481984 < .001980.

## Threshold-focused proposal
On old numerical near subspace Q, let A be the compressed regularized comparator, K the candidate metric after eliminating old carriers, R the residual cross factor on Q. Choosing j subtracts u u*, u=sqrt((8/35)/K_jj) R_j*. For invertible A-tau I, block congruence of [[A-tau I,u],[u*,1]] yields n_+(A-tau I-u u*)=n_+(A-tau I)+n_+(1-u*(A-tau I)^-1u)-1. Away from equality, a strict crossing is predicted iff u*(A-tau I)^-1u>1. Use tau=.99, prioritize predicted crossings, then positive-side weighted gain with weights 1/(1+((lambda-tau)/.05)^2) for lambda>tau, zero otherwise. Near-singular shifts disable crossing test. Update A,R,K by rank-one subtraction/Schur elimination. Floating-point proposal only. Match 32-carrier trials at 12 low and 20 high carriers from same shortlist. Threshold heuristic predicts 32 crossings (191 to 159 compressed eigenvalues above .99), not a certificate.

## Executed spectral acceptance
All four use BD ridge-absorbed error and exact-grid Gram/residual checks; all base/product errors retained. Subsets are principal restrictions of same enclosed checkpoint, retaining full operator error allowances.

| Trial | Carriers | Ridge | Absorbed | Even constraints incl constant | beta diagnostic |
|---|---:|---:|---:|---:|---:|
| BE inherited |832|1/256|4889|192|.9956839393|
| BH leverage |864|1/256|4905|176|.9949191750|
| BH threshold |864|1/256|4912|169|.9960264956|
| BH full |896|1/256|4920|161|.9948376864|
| BH full |896|1/512|4925|156|.9954433990|

Exact strongest beta=33750490965650593138084326586351832973829811433/33904982440892487352819727970703841835981209600 <997/1000. All four exact comparisons checked in bh_check.py. Threshold heuristic certifies seven more absorbed directions than leverage at matched composition and ridge; single test, no universal superiority or optimum claimed. Full batch at smaller ridge strongest. Even 5080-4925+1=156, inherit BG odd 151, total <=307 (36 fewer than BG). Physical margin >=3/2000. Upper bounds, not observed negative eigenvalue counts. Constant still excluded; no joint pole projection or new odd run.

## Reproduction
Use BD/BE even helpers and arrays and actual upstream derived matrices in separate writable BH directory. Verify BE manifest; link only read-only inputs. Run bh_select.py, bh_cross.py F/G0/G1/G2/G3, bh_metric.py new/cross/combine, bh_schur.py, bh_prebuild.py, bh_finish.py, bh_threshold.py, bh_baseline.py, then OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=256 python bh_ridge.py, bh_threshold_ridge.py, bh_leverage_ridge.py; repeat bh_ridge.py at RIDGE_DENOM=512. bh_check.py checks custody, preservation and recorded exact comparisons. Metric/check scripts use ../cr_int_t1d_be; select/cross accept BH_PARENT. Projected archive stores actual full effective metric, interval blocks, carrier IDs, incremental rows, near-space proposal and manifest. Extract with BH scripts and inherited helpers to rerun projected gate or regenerate subset inputs. Full upstream inverse needed only for full rebuild. Do not assume future local persistence.

## Next BI
Use threshold-focused selection for new batch on both parity checkpoints with actual marginal residual metric; test smaller ridge on odd BG. Keep matched composition comparisons. Joint pole-constrained inverse and constant coupling remain open; do not subtract dimensions without certificate. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.
