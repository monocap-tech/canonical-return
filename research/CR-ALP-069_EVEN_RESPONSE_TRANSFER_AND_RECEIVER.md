# CR-ALP-069 — frozen even residual inverse comparison PASS; original receiver open

DRAFT / NONPROMOTING. Same frozen even K1250 comparison witness from CR-ALP-056/058, same archived physical residual, same weighted CR-ALP-067 rational trial, same all seven signed prime shifts. Conditional on inherited classical source identities, CR-ALP-058 proxy and interval quadrature premises, and CR-ALP-068 pole-centered Gamma estimate. No Lean, full collar or RH claim.

## Custody

Recovery workflow 45621094974 job 139229123094 SUCCESS. It downloaded the original artifact 15023858000 from run 45620048012, not a recomputed solve. Original 982-line trial and diagonal files were committed/read back at 698c70e446c9197feff7588394e6b8c88c9b9679. Printed artifact SHA256s: trial 27670ee547ec2f5fc30fd3adb155378267e82df6f49bbaf90701c31f1fe4c15e; diagonal 38024951a1daf15f1a19c0298d793443458e9033e5b5620af3178a8378728126. Prior coefficient and residual logs and exact endpoint were retained at 9c8618cc7be833e68860fec507933e542d4154a4. Independent hash replay submitted, not yet passed at note creation.

## Exact algebraic transfer

Let x be the recovered exact projected 067 response trial, rho its archived physical residual proxy, D_old the old blended diagonal and D_new the full 068 diagonal, rounded down to 1e-9 per band. Define Delta=(D_new-D_old)x. Delta is a finite ordinary sine expansion through j=2048. Reconstruct x by the exact finite pole projection using saved outward pole coefficient intervals. Reconstructed X2 overlaps old native interval. For unprojected response residual y, compute E_new=E_old-<x,Delta>, y_new=y_old-Delta, norm2_new=norm2_old-2 Re<y_old,Delta>+||Delta||², pole_pair_new=pole_pair_old-<Delta,b>. The finite y_j update is exact algebraic subtraction; low sine entries unchanged. Use actual Q85 projection and Parseval for the entire tail. Early 068 band floors are nonmonotone; take the actual minimum over all remaining bands. For dmin=0.265358630, dtail=1.822938456, compute the weighted inverse residual bound. Retain 1e-13 translation transport and 058 source proxy L2 error. BigInt outward fixed point scale 10^30, integer roots upper, all inequalities exact. No new physical integration.

Results at 9b4b3935130e82966519214b57265d3fe154b021 (source and output read back):
- E_new upper 0.010156517621250860313072403446
- complete projected residual norm2 upper 0.086536930247193016254609837732
- infinite tail norm2 upper 0.058294025593250175541788017518
- weighted residual upper 0.053303622021144691967020696913
- proxy inverse upper 0.063460139642501891726798158073
- true frozen r inverse upper 0.063491442288544510816549782803
- old K lower 0.084343420037917007929466533262
- strict comparison margin 0.020851977749372497112916750459.

Thus <r_L,C^-1 r_L><K_L for the frozen even residual and old energy. Odd frozen residual already passed under L in 064. The even comparison is PASS_WITNESS in the limited sense of this inequality; do not conflate it with acceptance of the original corrected Schur form. Independent Python Fraction verifier reorganizes the norm update by projecting Delta first, checks all five recovered original artifact SHA256s, and asserts margin >0.02. Workflow submitted at e7dc5d81deb23b92cb85e1fbb57ecdeac7dca7bb; native independent replay pending.

## Critical original-form receiver obligation

CR-ALP-056 explicitly formed the original frozen correction X' against L=gamma I+P_Q: K_L=H-R*X'-X'*R+X'*L X', r_L=R-L X'. Write C=L+A, A=Q85(D_new-gamma I)Q85>=0. Completing the square against C instead gives K_C=K_L+X'*A X' and r_C=r_L-A X'. It is NOT generally valid to use K_L and r_L unchanged with C^-1 to assert original block positivity: the additional form A acts on the unshifted physical complementary variable, and shifting by X' creates mixed terms. In the old shifted coordinates the added term is (q-X'p)*A(q-X'p), not simply q*A q. The CR-ALP-069 comparison has a different receiver from the original completed-square gate; numerical PASS does not establish K_C-r_C*C^-1 r_C>0. To receive into original form, reconstruct the original X' and compute the correction energy and residual changes (including the complete pole bridge and tail), or derive another valid lower comparison for the whole shifted block. An alternative is an independently certified bound for the original unshifted Schur complement. These obligations must be resolved before claiming new original positive directions. Even the two scalar witness tests are not a full matrix certificate. No new 85-column positivity, index decrement, Phase IV or RH promotion. Prior 1048 positive dimensions and conditional cap 66 unchanged. Historical standalone archive reconciliation remains open.
