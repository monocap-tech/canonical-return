# CR-INT-T1D-BD — Full even quadratic metric and ridge-absorbed errors

**Status:** executed local directed-grid audits PASS. Conditional negative index <=707 (354 even + 353 odd), complement margin 3/2000. Original fixed window (-2,2), smooth double Mellin-pole-null source, all 24 signed prime powers. Inherited analytical source comparisons and ordinary binary64 exact-integer-dot arithmetic model. No Lean, remote CI, independent mathematical review, or new full frequency scan.

## Full even reconstruction
Adapt AZ's accelerated Gram kernel to positive even indices by adding the sum-frequency term; the constant j=0 remains separately excluded. All AQ/AS and AW source hashes pass after copying the initially omitted joint_final.log. Kernel comparison with inherited interval routines PASS; supplementary exact Python integer checks 23,699 entries (7,080 M; 2,000 C; 14,619 H). Deficits 5,080; surplus 12,619. Errors M 2711/67108864, C 4063/67108864, H 25/262144. First combined solve/base execution timed out after residual; separate base rerun completed. True residual 172258760381/140737488355328 < .001224; base error 706158969205005615/147573952589676412928 < .004786. Select 128 low +512 high carriers with AZ leverage heuristic, not globally optimal. Six cross blocks and combination PASS. Carrier solve residual 190446061/274877906944 < .000693.

## Exact full metric
c=8/35, J0=(399/50 I-A)/4, K=P*J0P/21+c(J0P)*(J0P). 49 physical cells, all 24 prime powers, 640 carriers; 409600 independent AX entry-enclosure overlap checks PASS. Max row interval width 53417/10^9. Full-space Gershgorin tests do not establish pilot percentages. Effective metric D=K-c Q_Y with Q_Y=GY+Y*G*-Y*HY. Errors eta_L=533521134964543/288230376151711744 < .001852, eta_D=103310254733405/252201579132747776 < .000410, eta_Q=462434685737425/288230376151711744 < .001605.

## Ridge absorbs matrix uncertainty
For arbitrary trial Z the exact bound is B <= B0_true-c(L*Z+Z*L-Z*DZ). Let Lhat,Dhat be enclosed grids with operator errors eta_L,eta_D. Construct F_r=B0_hat-c(Lhat*Z+Z*Lhat-Z*(Dhat+r I)Z), r>eta_D. For any x and y=Zx, the error relative to F_r, apart from base and product enclosures, is <= c[2 eta_L ||x|| ||y|| -(r-eta_D)||y||^2] <= c eta_L^2/(r-eta_D)||x||^2. Thus use the ridge already in the comparator to absorb L,D errors, instead of adding a separate global/projected trial norm penalty. Keep all base and directed product allowances. Valid for arbitrary Z, solver only proposes. Both parities satisfy r=1/128>eta_D. Exact scalar Young allowance included in bd_ridge.py, and inherited integer-grid Gram and residual test checks subspace.

## Executed spectral results
| Sector/trial | Absorbed | Constraints | beta diagnostic |
|---|---:|---:|---:|
| even standard 1/64 | 4705 | 376 | .9921486430138096 |
| even standard 1/128 | 4709 | 372 | .9969462026852349 |
| even ridge-absorbed 1/128 | 4727 | 354 | .9955528510767413 |
| odd ridge-absorbed 1/128 | 4790 | 353 | .9935559639738104 |

Even exact beta 19298223251037177230102589527869006698516606211685/19384428692223785688928813622131894949375490654208; odd exact beta 281201505932259672842950252897321084911826084206822379/283025331363892849495778610697886388712695178534584320. All four rational comparisons <997/1000. Odd input matches BC SHA256 9b0bd4864961c95f3ffaf3cda9b20470355942e02735626b6d2b780273e915ab. Odd inverse not rebuilt; new projected test executed. Even count 5080-4727+1=354; odd 5143-4790=353. Conditional index <=707, 278 fewer than BC; q0>=3/2000 ||f||² on certified complement. Upper bounds, not observed negative eigenvalue counts.

## Reproduction and custody
Local BD reproducibility archive contains even-adapted scripts (historical az_*, ba_*, bb_* filenames retained), source arrays, logs, manifests, bd_check.py. Separate even projected checkpoint contains actual bb_metric_input.npz, bb_metric_0.npz, az_carr.npy and hashes. With OPENBLAS_NUM_THREADS=4 and RIDGE_DENOM=128 run bd_ridge.py on even checkpoint. For odd, use BC's checkpoint in separate directory and bd_ridge_odd.py. Full upstream regeneration is possible with supplied scripts. Large intermediate inverse matrices are not included in projected archive. No guarantee of cross-session local persistence.

## Next — BE
Investigate targeted carrier expansion using residual directions and anisotropic ridge/error compensation on existing checkpoints, rather than just tuning scalar ridge. Treat constant even coordinate and original pole constraints explicitly. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.
