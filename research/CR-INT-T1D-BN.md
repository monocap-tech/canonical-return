# CR-INT-T1D-BN — Joint constrained signed-energy minimization and finite-head pilot

**Status:** local exact directed-grid audits PASS in both parities. Twenty complementary directions per parity (retain BM's four and add sixteen); exact joint energy-minimizing, four-moment-null lifts preserving every BK constraint. Conditional negative-index bound remains <=240 (121 even, 119 odd); BK coercivity margins 17/10000 and 53/10000 retained. A rectangular projected finite-head pilot through Neumann index 256 passes, but no whole-head or energy-normalized residual sign certificate. Conditional on inherited fixed-window (-2,2), double Mellin-pole-null form identities, AF/BL complete Neumann action, AM arithmetic norm, BK complement and form/core transfer, BM tail theorem. No new frequency scan, Lean, remote CI or independent mathematical review.

## 1. Joint optimization

Let V be BL's selected d=121/119 pole-null Neumann lifts. Extend V_n from four to twenty modes using a numerical greedy log-determinant proposal on the normalized moment-response matrix, candidate unused same-parity Neumann indices 200+p through 5199. The proposal is not a certificate; selection need not be globally optimal. All four old modes are preserved. For C=(F*V)^-1 F*V_n, U=V_n-V C, F*U=0 exactly. Recompute directed BL cross overlaps; old block preserved. Evaluate the original signed energy Q on [V,V_n] using BL's high-precision complete Gamma divided differences and exponential tail, all 24 signed prime powers, grid 10^-110. Original 125/123 square blocks preserved exactly. Direct interval LDL again inconclusive, but integer congruence positivity PASS for 141/139-dimensional extensions (.999I < Z Q Z* < 1.001I). Widths <3.46e-73 /3.49e-73 respectively. No inference of negativity from the failed LDL.

Define J=U*L V, D_U=U*L U, and BM normalized moment rows M_o on V, K on U (normalization 1,1/4096,1/4096^2,1/4096^3). Exact joint response X minimizes signed energy of W=V+UX subject to K X=-M_o. Saddle system [D_U K*; K 0][X;Lambda]=-[J;M_o]. The audit scales both K blocks and moment RHS by 2; multiplier coordinates rescaled accordingly. D_U>0 by finite positivity and independent U columns, so a solution is the unique constrained minimizer for each exceptional coordinate. Since four-mode BM trial is feasible in enlarged space, exact minimum energy cannot exceed its energy (in Loewner order). This is only an upper approximation to full physical Schur.

## 2. Exact residual and coefficient bound

Use S=10^12 rational-grid centers, Python arbitrary-precision integer products, and Frobenius upper bounds obtained by exact integer square-root rounding. Let P0 be central scaled saddle matrix and P exact one; let eta_P bound its norm error. Likewise B0 is RHS center with error eta_B. Numerical inverse Y and solution X0 are proposals. Verify rho=||Y P0-I||_F+||Y||_F eta_P <1/2, and ||X_exact-X0||_F <= ||Y||_F (||P0 X0-B0||_F+eta_P||X0||_F+eta_B)/(1-rho). Include all grid rounding allowances. Exact correction coefficients in the pole-null basis are [I-C X;X] (X here the first twenty rows of the saddle solution); reference mode coefficient from the exact pole-null ratios. The stored coefficient_error is the Frobenius error for these non-reference coefficients. BM tail coefficient 2->1 norm includes the reference coefficient. Results (diagnostics, full exact fractions in actual NPZ):

| Quantity | Even | Odd |
|---|---:|---:|
| C error | .000189579 | .000280718 |
| J error | .011569037 | .018696011 |
| D_U error | .059881774 | .092843968 |
| K error | .002555489 | .003592019 |
| rho | .066545135 | .103665425 |
| saddle solution error | .383027899 | .539171408 |
| correction Frobenius bound | 3.541062 | 2.860276 |
| coefficient error | 1.333039 | 1.863989 |
| full coefficient l1 operator bound | 155.651942 | 132.400149 |

All accepted comparisons use exact rational residual bounds, not numerical solver accuracy. The large coefficient enclosure errors remain important for near-null energy normalization.

## 3. Infinite signed-action tail

BM's exact M0,H0,M1,H1 cancellations hold identically for the exact minimizer, not just the proposal. All retained frequencies below W=4096, so reuse BM tau_N=(74 W^4 A)^2 (4/3)^10/(18 N^9), with the new certified A values. Relative to BM's stored A, improvement in squared prefactor is about 528854 even and 16329 odd. All exact checks PASS:

| N | even tau_N upper | odd tau_N upper |
|---:|---:|---:|
| 2^18 | 2e-12 | 2e-12 |
| 2^20 | 7e-18 | 5e-18 |
| 2^24 | 1e-28 | 1e-28 |

No enumeration to these cutoffs. This is unprojected operator-action tail in original BL coordinates, not a bound on the tail of P_C L W or on energy-normalized residual. BM projection-safe weighted inverse remains available.

## 4. Rectangular finite-head pilot

Evaluate signed rectangular mixed action on orthonormal Neumann indices m<=256 of the correct parity, with full diagonal where m=n. For off-diagonal nonexponential part use 2(d(w_m)-d(w_n))/(w_m^2-w_n^2), d as in BM; subtract reference mode contribution to make columns exactly pole-null. Include signed exponential Gamma k>=1 terms and directed geometric tail; k=0 vanishes on pole-null columns. Divide rows by sqrt(physical Neumann mass). Directed grid 10^-12 produces R0 and error e_R. With joint coefficient center B0 and e_B, central head Y0=R0 B0, ||Y-Y0||_F <= ||R0||_F e_B+e_R(||B0||_F+e_B). Use independently constructed physical F*head overlaps and inherited G=F*F enclosure. For arbitrary projection trial T, ||P_C Y0||_F <= ||Y0-F T||_F, and square = ||Y0||_F^2-2 Re tr(T* F*Y0)+tr(T* G T). Numerical T approximates G^-1 F*Y0; exact integer products and directed overlap/Gram errors certify upper bound independently of optimality. No need to project the tail or assume commutation. Add ||Y-Y0||_F after projection. The trial projects only onto F, not additionally onto g, so is conservative.

| Quantity | Even | Odd |
|---|---:|---:|
| Pilot rows |129|128|
| rectangular action error |1.35e-10|1.34e-10|
| physical overlap error |5.91833e-6|1.101552e-5|
| central head norm |23.140207|26.024207|
| central trial projection bound |17.533824|17.019260|
| lift uncertainty allowance |33.947331|50.221936|
| final projected head upper |51.481155|67.241195|
| unprojected upper |57.087537|76.246143|

All entries in unnormalized exceptional coordinate Frobenius metric. Cross-check rectangular action after pole-null row congruence against independently built square original signed energy: 8037 even,8340 odd overlap checks PASS. The central residual reduction is genuine in the certified upper comparison, but the coarse lift uncertainty dominates. The cutoff 256 is only a pilot, below BM's tail validity threshold; it cannot be combined with the displayed large-cutoff bounds. No full head, original energy-normalized residual or Schur sign established. Do not subtract any exceptional count.

## 5. Reproduction and next BO

All 10 BM, 45 BL and 22 BK manifest entries verify. Archive includes actual new indices, interval matrices, rational preconditioners, saddle proposals and error enclosures, rectangular action and projection trials, scripts/logs and compact upstream manifest files/data. Run from cr_int_t1d_bn with sibling BM/BL/BK directories. Set OPENBLAS_NUM_THREADS=4; run bn_select.py; for p=0,1 run bn_lifts.py p, bn_energy.py p --extend, bn_precondition.py p --extend, bn_joint.py p, bn_head.py p 256; then bn_check.py. Selection is exploratory; stored indices can bypass it. Full upstream provenance requires the earlier chain; this archive supports BN-stage replay from inherited actual compact inputs. No independent formal verification.

**Next BO:** tighten physical constraint overlaps and joint lift coefficient enclosures with arbitrary-precision directed grids, preferably compute residuals in the signed-energy-preconditioned metric directly instead of multiplying the existing absolute error by a near-null preconditioner. Extend rectangular head evaluation or use structured analytic sums for the intermediate head; do not combine index-256 pilot with N>=2W/(3/4) tail theorem. Continue optimizing complement response beyond moment constraints. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.