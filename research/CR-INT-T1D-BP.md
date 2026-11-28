# CR-INT-T1D-BP — Full signed action and first direct Schur-positive subspace

Status: executed local audits PASS for complete finite head and projection, and for a three-dimensional odd Schur-positive subspace. Conditional negative index <=237 =121 even +116 odd. No unconditional positivity. Inherits BO/BK/BN analytical source comparisons, form/core hypotheses, and ordinary binary64 exact-integer-dot model for the projection products. No Lean, remote CI, independent mathematical review or new frequency scan.

## Complete rectangular action
For the BO energy-normalized lift coefficient matrix B, write the action as R+XT+E, where R is the complete signed nonexponential divided difference (with the exact Mellin pole reference subtracted), X and T factor exponential Gamma terms k=1,...,39 and E is bounded by the inherited geometric remainder. Contract R B+X(TB) on directed 10^70 grids, with rounded coefficient grid 10^65 and final output grid 10^12. Account for all coefficient, action, factor and output errors via Frobenius bounds. Reuse BO's checked head through 256. Tiles of 256 rows are checkpointed; phase recurrence reseeded every 32 steps. For high frequencies m>256 off the selected trial mode set, direct Euler–Maclaurin at z=1/4+i*w_m/2 avoids 96 recurrence steps. With x=Im z>100, |z+t|>=max(x,t), so the periodic Bernoulli integral remainder is at most |B_64|(65/64)x^-64. Directed shifted-EM overlaps and direct phase checks PASS. Full head at common cutoff 16384: 8193 even, 8192 odd rows. Certified action errors <9.92e-7 and <9.75e-7; exact rational bounds are stored in the head checkpoints. Frobenius upper 22.731689/22.484219. 1536 direct phase overlap checks each. Independent old-style directed complete signed row calculation on six new rows per parity agrees with full head enclosure (726/714 normalized entries). All 24 signed prime powers retained. No need to enumerate beyond cutoff: BO tail squared norm <1.32e-7/1.14e-7. Coverage complete, not yet sufficient full Schur sign.

## Full physical projection
Let Y be the actual orthonormal finite action head, Y0 its rational center. F is BK's pole-null exceptional synthesis, G=F*F, and A=F*|head. Enclose A on 10^12 grid, using inherited directed raw overlaps and exact pole subtraction. Exact radix-2^16 integer dot products via ordinary binary64 GEMMs (absolute dot bound checked <2^53); short products use Python integers. Use arbitrary trial V for the projection, proposed by a numerical solve. For all z, ||P_C Y0 z|| <= ||Y0 z - F V z||; not necessary to subtract the pole component of Y separately. For T0=A0 Y0, center P0=Y0*Y0-T0*V-V*T0+V*G0 V, with directed floor errors. Full overlap Frobenius errors <5.484e-5/1.089e-4, global Gram error <.014877/.034080; projected trace upper <427.969/388.105 before adding action error. Uniform scalar test absorbs no dimension. These are residual bounds, not original energy eigenvalues.

## Restricted Schur certificate
Let Q be a rational-grid proposal with ||Q*Q-I||<=delta<1/1000. Restrict the physical projection uncertainties before estimating norms:

  E_Q=2||Y0 Q||_F eta_A ||VQ||_F+eta_G||VQ||_F^2+rounding allowance,
  h_Q=(||Q*P0 Q||_row+E_Q)/(1-delta).

For every z in range Q, ||P_C Y0 z||<=sqrt(h_Q)||z||. The full dual energy residual is bounded by (sqrt(h_Q)+eta_Y+sqrt(tau_N))/sqrt(mu), where mu is BK's parity complement coercivity and tau_N BO's squared action tail. Projection need not commute with the head/tail projection; the triangle inequality is applied before taking the dual norm. BO's normalized lift energy is >=(1-eps)I. Thus strict rational test (sqrt(h_Q)+eta_Y+sqrt(tau_N))^2 < mu(1-eps) gives positivity of the exact original Schur form on range Q, conditional on BK's Schur reduction. Square roots rounded upward, Gram row norms and restricted error products computed exactly. Proposal eigensolver is not trusted.

| Trial | Full squared residual bound (rounded up) | Threshold (rounded down) | Result |
|---|---:|---:|---|
| even 1 | .024766 | .0016999984 | FAIL |
| odd 1 | .000820625 | .0052998847 | PASS |
| odd 2 | .002494779 | .0052998847 | PASS |
| odd 3 | .003489316 | .0052998847 | PASS |
| odd 4 | .011506615 | .0052998847 | FAIL |

Odd 3-dimensional actual rational proposal saved; restricted projection error <.001793, central projected norm-squared bound <.001658. Schur lower bound on this subspace exceeds .34 times normalized coordinate norm squared (from displayed conservative bounds). Inertia reduction: the original complement is coercive and the exact Schur matrix has a 3-dimensional positive subspace in its 119-dimensional odd sector, hence n_- <=116 odd and <=121 even, total <=237. Original BK margins on original complements unchanged; no new global coercivity margin claimed. Failed tests do not establish negative energy or optimality. This is a direct signed Schur reduction, not an additional carrier expansion.

## Reproduction and custody
Run in sibling directory next to BO, BK, BN inputs: OPENBLAS_NUM_THREADS=2 python bp_head.py p 16384; bp_project.py p; bp_restricted.py p; bp_check.py (p=0,1 for first three). 90 inherited digests checked (34 BO,22 BK,34 BN), full head coverage and supplementary old-style action overlaps PASS. Archive actual head matrices, projected Gram centers, trial coefficients, exact restricted subspace, scripts and logs. Separate upstream BO/BK/BN archives required for full reproduction. Tiles/precision/selection are experimental choices; no global optimum claimed. Initial invocation interrupted, actual complete runs PASS. Uniform projected Gram error and even/odd fourth tests fail; no further reduction claimed.

## Next BQ
Use the complete action to build residual-driven complementary correctors in normalized coordinates, prioritizing the unresolved Schur directions. Use a stronger frequency-dependent inverse or variational residual norm where justified; preserve exact constraint and tail moment conditions. Consider removing the pole component from the head projection as an additional valid improvement. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain unresolved.