# CR-INT-T1D-BO — Energy-relative precision and structured infinite tail

Status: local directed integer audits PASS for both parities. Conditional negative-index upper bound remains 240 (121 even, 119 odd). BK complementary margins 17/10000 and 53/10000 retained. No full residual or Schur sign established. All results conditional on inherited source identities, smooth-core transfer, and upstream numerical certificates. No Lean, remote CI or independent review.

## High precision overlap and constrained solve
Use exact sine values at integer multiples of pi/4 and directed sqrt/pi on the 10^75 grid; multiply the rational exceptional basis using arbitrary-precision integer products. The 17061 even and 16541 odd cross entries are compatible with old interval bounds. Cross error <1.4e-71/2.2e-71. Compute the same BN constrained minimizers on 10^80 grid, using five Newton inverse iterations from numerical proposals and exact residual acceptance. H inverse residual <7e-71/1.2e-70, KKT residual <8.5e-68/1.6e-67, corrected coefficient error <1.6e-66/2.5e-66. Four exact BM moment constraints and F*W=F*V preserved. Unnormalized coefficient l1 bounds <140/110.

## Actual corrected energy normalization
Let B denote exact corrected lift coefficients in the pole-null Neumann basis, A original signed finite energy and B0 rational center with error eB. Enclose A0=B0* Acenter B0 with error eta <= rounding + ||Acenter||(2||B0||eB+eB^2)+||Arad||(||B0||+eB)^2. All norm bounds via exact integer Frobenius products. Fixed-point triangular LDL proposes invertible Z on 10^110 grid; check Z A0 Z* relative to identity by exact row sums, add eta ||Z||_F^2. Both actual corrected energy congruences pass: relative error <9.12e-7 even, <2.18e-5 odd; normalizer norms <1.804e28/7.050e28, normalized lift coefficient errors <2.86e-38/1.73e-37. Write W_hat=W Z*. No claim that the exact infinite Schur matrix is positive.

## Normalized finite-head pilot
Rebuild full signed rectangular action at index <=256 on 10^75 grid using inherited directed Gamma and prime evaluator; reuse BN's physical F*head overlaps and BK's physical Gram enclosure. Exact integer products produce the normalized central head. For any trial T, ||P_C Y||_F <= ||Y-FT||_F, checked by exact quadratic expansion with directed errors. No projection of the tail. Even/odd action errors 9.04e-73/9.29e-73; propagated lift/action errors <7.27e-37/4.65e-36. Central head Frobenius norms 13.7814613/14.7080122, projected bounds <12.408953/12.625282. This is not small enough for scalar Schur acceptance. Supplementary rectangular/square enclosure overlap checks 8037/8340 PASS. Initial attempts overflowed when converting large fixed-point integers before dividing; corrected by conversion of the exact rational ratios. Index-1024 extension timed out; no result accepted for it.

## Structured normalized tail
BM gives exact action pairing at exterior frequency x=w_m. Let W0=4096, r=40, M_l and H_l as in BM; M0,H0,M1,H1 vanish exactly. For normalized lifts, use actual row norm enclosures m_l>=||M_l(W_hat)||/W0^(2l), h_l>=||H_l(W_hat)||/W0^(2l+1), l=2,...,39. Compute the entire signed row before taking its Euclidean norm, including pole reference and exponential terms. For exponential moment terms k>=40, ratio bounded by exp(-8)((a40+2)/a40)^(2r)<1/2; bound omitted normalized moments by twice first term. Exact check sum_{k>=1} exp(-4a_k)(a_k/W0)^(2r)<1. With A the full normalized coefficient 2->1 bound, x>=2W0, finite geometric division of denominators gives

    |<p_m,L W_hat z>| <= ||z|| [sum_{l=2}^{r-1}(36 m_l W0^(2l)/x^(2l+1) + 2 h_l W0^(2l+1)/x^(2l+2)) + 74 W0^(2r) A/x^(2r+1)].

Here |d(x)|<=18x and |d(w)|<=18W0, x^2-w^2>=3x^2/4. Nonexp remainder <=48 W0^(2r) A/x^(2r+1)+48 W0^(2r+1) A/x^(2r+2); exp remainder <=2 A W0^(2r)/x^(2r+2); sum <=74 W0^(2r) A/x^(2r+1). No convergent expansion in unbounded a_k assumed. Evaluate bracket B_N at x=3N/4. All terms decay at least (N/m)^5, so squared L2 tail (mass 2) <= N B_N^2/18. All 38 pairs of signed moment rows enclosed separately. No new moment constraints beyond first four. Example normalized M2,H2 norms <.282/.588 even and <.264/.525 odd; the huge individual coefficient norms cancel in these rows. Exact bounds:

| N | Even squared tail | Odd squared tail |
|---:|---:|---:|
| 16384 | 1.32e-7 | 1.14e-7 |
| 65536 | 4e-13 | 4e-13 |
| 262144 | 1.4e-18 | 1.3e-18 |
| 1048576 | 5.3e-24 | 4.6e-24 |

Contrast crude BM-style normalized bound at 16384, approx 2.394e53/3.700e54 using the same A. These are upper bounds, not measured tail norms. With BK scalar inverse margins, normalized dual squared tail at 16384 <7.77e-5/2.16e-5. Finite head through 256 cannot be combined with tail beyond 16384: intermediate frequencies unaccounted. Also even full head norm would need sufficient complement response and tighter inverse comparison. Conditional index stays <=240.

## Reproduction and next BP
111 inherited SHA256 checks PASS (34 BN,10 BM,45 BL,22 BK). Run scripts in BO directory with sibling BN/BM/BL/BK inputs. OPENBLAS_NUM_THREADS=4 python bo_lifts.py p; bo_joint.py p; bo_normalize.py p; bo_head.py p 256; bo_tail.py p; bo_check.py. p=0 even,1 odd. bo_head_extend.py prepared and attempted at 1024, timed out. Retain all new actual normalized coefficient grids, moment bounds, scripts and logs, with SHA256 manifest. Parent BN archive carries upstream dependencies; full provenance needs the earlier chain. No new trial modes, no new frequency scan. Choices of grid, expansion order and cutoffs were made for this experiment. BP: accelerate rectangular action by separable exponential factors, tiled nonexponential signed divided differences and exact products; bridge 256 to 16384; then target additional complement responses in normalized energy coordinates, since current head upper bound is not enough. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.