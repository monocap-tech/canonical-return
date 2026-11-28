# CR-INT-T1D-BA — Carrier-space joint surplus reoptimization and projected error

**Status:** conditional local directed-grid certificate, fixed window I=(-2,2), original smooth double Mellin-pole-null source. All 24 signed prime powers retained; Gamma image correction still unused. No full positivity, Lean, remote CI, or independent mathematical review.

## 1. Carrier-space reduction

Reopen AZ and AY. With c=8/35, H=D0+R*R, C=T*R, E_V=HV-C*, B_V=M-CV-V*C*+V*HV, F=P*J0T, G=P*J0R, L=F-GV, D_P=P*J0P+epsilon I, the inherited trial upper bound is B_V-c(L*Z+Z*L-Z*D_P Z). Set V'=V-cYZ for arbitrary Y. Expansion, including the change L'=L+cGYZ in the arithmetic cross term, yields

    B <= B_V-c(L_Y*Z+Z*L_Y-Z*D_Y Z),
    L_Y=L+Y*E_V=F-Y*C*+(HY-G*)*V,
    D_Y=D_P-c Q_Y,
    Q_Y=GY+Y*G*-Y*HY.

No numerical solver correctness is assumed. Choosing Y=H^-1G* maximizes Q_Y in Loewner order by completion of squares; if E_V is nonzero it need not optimize the whole trial for a fixed Z. The family contains the old trial at Y=0; arbitrary choices of Y,Z need not improve it. For the optimal metric choice, Q_Y=GH^-1G*. This requires 640 right-hand sides instead of a fresh 5143-column inverse. All terms are finite on the retained comparator; no exact full arithmetic inverse is claimed.

## 2. Executed odd checkpoint

All 14 entries of az_derived_sha256.txt checked against the local archived matrices. The first combined attempt ba_schur_interrupted.py timed out after the Y residual check; no spectral acceptance from that run. The split ba_schur.py and ba_metric.py both completed. Y rounded on 2^21 grid; original H,C on 2^26, arithmetic F,G,D on 2^25. Exact integer products use AR ordinary IEEE binary64 absolute-dot range checks and radix splitting where required. The residual E=HY-G* is enclosed on 2^26 grid, with error <= eta_H ||Y||+eta_G+sqrt(nr)/2^26. The exact bound is

    ||Y|| <= 6187089/2097152,
    ||E|| <= 33944710515/35184372088832 < 965/1000000.

For F-Y*C*+E*V the error is eta_F+||Y|| eta_C+eta_E ||V||+2 ceil(sqrt(rm))/2^26. For Q_Y=sym(GY)-sym(Y*E) the error is eta_G ||Y||+||Y|| eta_E+4r/2^26. Sym denotes the arithmetic average of a matrix and its adjoint, with explicit integer rounding. The effective D is stored without ridge, to be added exactly for each new proposal. Enclosures are

    eta_L = 111872791660381735/36893488147419103232 < .003033,
    eta_Q = 171903217983166819/73786976294838206464 < .002330,
    eta_D = 187306001498773859/322818021289917153280 < .000581.

The same inherited AZ base error eta_base and 640 carriers are used. No even reoptimization, new carrier set or new frequency scan.

## 3. Projected spectral verification

The uniform-error postprocessor ba_spectrum.py at ridge 1/16 passed with 4551 absorbed dimensions. For better use of the error budget, ba_projected.py keeps the fixed matrix error eta_fixed (base plus directed trial-product rounding) separate from the L,D input errors. For the fixed-point proposed basis Q, let delta>=||Q*Q-I||, eps>=||Mgrid Q-Q Lambda||, ell>=||Lambda|| and zq>=||ZQ||. The upper Rayleigh quotient is bounded by

    beta = ((1+delta)(ell+eps+eta_fixed)
            + c(2 eta_L (1+delta) zq + eta_D zq^2))/(1-delta).

Here ||Q||<=1+delta and ||Qx||^2 >= (1-delta)||x||^2. The fixed product errors were bounded globally, but L,D errors are paid only in the proposed subspace. ZQ is calculated by checked exact integer products followed by directed floor on the 2^21 grid, with ceil(sqrt(r dim))/2^21 Frobenius error added; its norm is bounded via a small exact Gram. A broad nested eigensubspace is proposed numerically, and Gram/residual bounds for it also bound every initial-column restriction. Each trial uses exact rational strict beta<997/1000. No numerical eigenvalue alone is accepted.

| Method | Ridge | Absorbed | Odd constraints | beta diagnostic |
|---|---:|---:|---:|---:|
| Uniform | 1/16 | 4551 | 592 | .995674656 |
| Projected | 1/16 | 4555 | 588 | .995349512 |
| Projected | 1/32 | 4562 | 581 | .996712583 |
| Projected | 1/64 | 4563 | 580 | .995938599 |

Strongest exact beta = 415075479277690024146568166971736861883999/416768141817570046685465257727520617267200 < 997/1000. Restricted norm <=10775507/2097152; global norm diagnostic 13.637875. The ridge values and eigenspace proposal were exploratory choices, not optimized globally.

## 4. Combined conditional result

Odd codimension <=5143-4563=580. Combine with AY's separately certified even codimension 608, inherited unchanged. Thus n_-(q0;D^0)<=1188 (14 fewer than AZ's 1202), with q0[f]>=3/2000 ||f||^2 on the certified complement. These are codimension bounds, not actual negative eigenvalue counts. The arithmetic and analytical assumptions are inherited from AM/AQ/AS/AV/AW/AX/AY/AZ. No full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV or RH proof.

## 5. Reproduction

The BA archive includes the AZ reproducibility scripts and logs, not the large derived matrices; regenerate AZ per its note. In the regenerated AZ directory, verify original arrays and run AZ's check, then run:

```
export OPENBLAS_NUM_THREADS=4
python ba_schur.py
python ba_metric.py
RIDGE_DENOM=16 python ba_spectrum.py
for r in 16 32 64; do RIDGE_DENOM=$r python ba_projected.py; done
python ba_check.py
```

The checkpoint hashes record the executed local arrays; other floating proposals may differ and must be reverified. The interrupted script/log are not acceptance evidence. BA reuses AZ's actual enclosed inputs (not a fresh scan or re-audit of the entire original chain).

## 6. Next — BB: remove the worst-case arithmetic metric

Current trial uses z=B^(1/2)PZ/21 and B+B^2<=21B. Instead keep the full metric K_P=P*(B+B^2)P/21^2. With c=8/35, scaled cross cG, the joint optimized carrier metric is K_P-c^2 GH^-1G* (in the actual quadratic cost, before factoring out c). K_P can be obtained from Gram pairings of B applied to the carriers, retaining the signed prime shifts, not from an unverified norm floor. Alternatively bound the uncertainty anisotropically by 2 eta_L ||x|| ||Zx|| <= a||x||^2+(eta_L^2/a)||Zx||^2, producing a verified matrix correction rather than a global norm budget. These proposed next computations have not been executed. Full fixed-window/global positivity remains open.