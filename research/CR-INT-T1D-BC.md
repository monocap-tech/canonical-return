# CR-INT-T1D-BC — Full quadratic arithmetic metric and recovered checkpoint

**Status:** executed full odd-sector spectral audit PASS. Conditional negative-index upper bound 985 (608 inherited even + 377 odd), complement margin 3/2000. Original fixed window (-2,2), smooth double Mellin-pole-null source, all 24 signed prime powers. Inherited analytical and exact ordinary binary64 integer-dot arithmetic assumptions; no Lean, remote CI or independent mathematical review.

## Checkpoint recovery

Regenerated AZ using the BB archive: az_prepare.py, az_spotcheck.py, az_solve.py, az_base.py, az_cross_checkpoint.py select and F,D,G0,G1,G2,G3, then az_combine.py. All 14 AZ digests match the historical az_derived_sha256.txt exactly; 23,954 supplementary integer spot checks PASS. Base residual 238283983419/140737488355328 < .001694. ba_schur.py regenerated the original carrier solve byte-for-byte; residual 33944710515/35184372088832 < .000965. Duplicated ba_metric.py timed out, so its historical effective-metric digest was not reproduced; bb_build.py directly substitutes the stronger metric and passed. AQ/AS and AW manifests passed.

## Exact full metric

For c=8/35, J0=(399/50 I-A)/4, K=P*J0P/21+c(J0P)*(J0P). bb_metric.py az_carr.npy passed on all 640 actual odd carriers across 49 cells, with 409600 independent AX entry-enclosure overlap checks and maximum row interval width 52237/10^9. Full-space Gershgorin relative tests did not certify the pilot percentages; no such extrapolation is made. Effective normalized metric K-c Q_Y, where Q_Y=GY+Y*G*-Y*HY. bb_build.py passed with eta_L=111872791660381735/36893488147419103232 < .003033, eta_D=185473253054224739/322818021289917153280 < .000575.

## Executed projected spectral audit

All three bb_projected.py runs passed the strict rational beta<997/1000 test, using BA's projected error and checked exact-grid Gram/residual products:

| Ridge | Absorbed | Odd constraints | beta diagnostic |
|---|---:|---:|---:|
| 1/32 | 4752 | 391 | .9909914393 |
| 1/64 | 4763 | 380 | .9918926013 |
| 1/128 | 4766 | 377 | .9964694900 |

Strongest exact beta=5612070728173517212807056830571177774291/5631954399664197343846725151056147251200. Odd improves from 580 to 377; even unchanged at 608; total conditional index <=985, physical margin >=3/2000. These are upper bounds, not actual negative eigenvalue counts. Ridges are exploratory, not globally optimal.

## Reproduction and custody

Local reproducibility ZIP includes inherited arrays, scripts, manifests and logs, including bc_check.py which checks local custody and exact recorded rational acceptance. Separate projected checkpoint ZIP includes actual bb_metric_input.npz (SHA256 9b0bd4864961c95f3ffaf3cda9b20470355942e02735626b6d2b780273e915ab), bb_metric_1.npz and original az_carr.npy. Unpack both archives and run OPENBLAS_NUM_THREADS=4 RIDGE_DENOM=128 python bb_projected.py to rerun projected gate. Full upstream regeneration is needed for complete provenance checking, not for projected rerun. All large regenerated matrices also retained in local working directory, whose cross-session persistence is not guaranteed. No full frequency rescan.

## Next — BD

Mirror full exact metric in even sector; generalize AZ's accelerated Gram kernel and record actual carrier and effective-metric checkpoints. The odd projected checkpoint can be reused for an anisotropic matrix Young-bound refinement. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.
