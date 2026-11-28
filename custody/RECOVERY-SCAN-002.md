# Recovery scan 002 — original fixed acceptance recovered

BP/BS/BT/BU uploads received; all ZIP CRC tests passed. SHA-256:
- BP 8,887,439 bytes: e82efaa2bd8360221d574e43d0fa5eb55f7de281f180d83eca6f1adc951624be
- BS 8,112,156 bytes: 892a7c143a79cf87c608dc3bf3eb89f30b8f3f7f33c2864a83d25ce24559ff8c
- BT 4,884,132 bytes: a67ec45b190f9de3a69aca2abcaa933229efaa2371f1672a4f9f40848c7f2aea
- BU 35,864,643 bytes: 55ab94bd18a384b019fd66d0de0989076b0a1a726556469db1418fab93800bb2

159 additional exact hash matches; total 218/578, 360 missing, zero mismatches. All five immediate blockers restored. All 166 CA/CB local manifest checks and five independent terminal comparisons, seven matched arrays PASS.

Original fixed rational scripts executed using recovered actual inputs, OPENBLAS_NUM_THREADS=4:
- cr_int_t1d_cb/cb_matrix.py --fixed: parameters 1/50,4,1/175; exact lower diagnostic 0.04506999999971855; fixed PASS (exact saved rational equality).
- cr_int_t1d_ca/bx_restrict.py --fixed: normalized lower diagnostics 0.9229629864206749, 0.5498874526488997, 0.2854024086714804; exited 0.
- cr_int_t1d_ca_odd/bv_joint_trial.py odd --fixed: direct 11D 0.1530776987278319, optimized 11D 0.12306684514700258, both fixed PASS. Failed 12D CA trial reproduced, not a negative Weil direction.

Actual local log SHA-256: CB add136b3dece344250fe02c9055dec04643000a8691d8cd94dadacb498e28b90; CA even 5bc4ffa4baaeec92690053c044ecd72ad38bcae7f40f0869915388fedfd53daf; CA odd 114a8e0a6fe2dfd9f934a26957d81e2d44fbd41d608b1deded862cffebad1e8d. Full local recovery report SHA-256 3e69749b5269e25f802ce73eb39158d4051dc2e9ee556e84ae20385666df73d8; Git status TSV is a reduced derivative without local origin paths.

No source/action regeneration, no remote CI, no CC acceptance. Archive bytes remain outside Git. Remaining recorded paths by stage: BI 2, BK 22, BN 33, BO 33, BQ 25, BR 25, BV 60, BW 31, BX 28, BY2 55, BZ 46. Conditional index <=225 unchanged. Original fixed check recovery is complete for CA/CB; full transitive source custody remains open.
