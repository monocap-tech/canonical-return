# CR-INT-T1D-AZ — Complete joint odd-sector arithmetic certificate

**Status:** local arithmetic certificate, conditional on the original source identities and the AM/AQ/AS/AV/AW/AX/AY analytical and data chain. Fixed window I=(-2,2), smooth double Mellin-pole-null source, all 24 original prime powers retained with original signs. No Lean, remote CI, independent mathematical review, fixed-window positivity or RH proof.

## 1. Full odd inverse, without a new frequency scan

The joint odd deficit coordinates number 5143 (=452+4691); the positive surplus coordinates number 12811 (=60+12751). The inherited input array hashes and low-frequency census passed again. AY's even result (608 constraints) is inherited, not regenerated. The objective is to replace AX's old odd count 724 by a single joint low/high certificate, rather than separate mixed-block tests.

Use AY's notation: t=1/2, D=I+R D0^-1 R*, D0=diag(1-alpha), H=D0+R*R, C=T*R, M=T*T, B=(6/5)(399/50 I-A), J0=B/(24/5). On the original physical source, q0[f] >= t(<(D+B)f,f>-||T*f||^2), with the nonnegative Gamma image term still unused. For any trial V, W=T-RV and E_V=HV-C*, the surplus trial upper comparator is B_V=M-CV-V*C*+V*HV, with gap to the exact retained-surplus inverse E_V*H^-1 E_V. From I/24<=H<=I the gap is at most 24||E_V||^2 I. This controls the retained surplus inverse only, not the full arithmetic inverse.

## 2. Accelerated exact grid and residual

The direct run of the inherited AY odd builder timed out during H construction; its partial log `az_odd.log` is not acceptance evidence. `az_fast_gram.py` implements the same directed integer-grid endpoints specialized to odd parity, using compiled Numba int64 loops. For each pair j,k the unweighted physical Gram is

    G_jk = (2/pi)(sin(pi(j-k)/4)/(j-k) - sin(pi(j+k)/4)/(j+k)),

with the difference term on the diagonal equal to 1/2. Since j,k are positive odd, the difference and sum terms have at most one nonzero off-diagonal sine; each nonzero denominator has absolute value at least 2. Thus the unweighted interval endpoints have magnitude <=10^9 on the original 10^9 grid, including the diagonal. The integer endpoint constants enclose 2/pi by the inherited Machin check. Weights are enclosed by integer square roots on the 10^8 grid. `checked_block` verifies bounds on integer products before the compiled kernel, which computes lower/upper endpoints on the 2^26 grid, checks every midpoint fits signed int32 and every entry error fits uint8, and returns the error arrays for exact row/column accumulation. Symmetry is checked on the full square grids. The initial small blocks also agree exactly with the inherited AY implementation, and `az_spotcheck.py` independently checks 2000 seeded entries of each grid and all M/H diagonal entries with unbounded Python integers (diagnostic spot check, not a replacement for the full kernel).

| Enclosure | Exact bound |
|---|---|
| M | 1291/16777216 |
| C | 8119/67108864 |
| H | 3203/16777216 |

All square grids and the rectangular cross grid were actually computed. `az_solve.py` proposes V by Cholesky; no correctness of that solver is assumed. `az_base.py` verifies integer grid dots under the inherited ordinary binary64 exact integer-dot model and uses the exact Gram V*V to tighten the operator-norm bound for V. It obtains ||V||<=3856253/1048576, instead of a loose rectangular row/column estimate. The actual full residual enclosure is

    ||HV-C*|| <= 238283983419/140737488355328 < 1694/10^6.

The exact retained-surplus inverse gap is <=24 times the squared bound, <7/100000. The symmetric grid approximation to B_V is built from M-sym(CV); the omitted sym(V*E_V) and all grid/product errors are included in its certified norm error

    eta_base = 1007208218832374735/147573952589676412928 < 683/100000.

This is a fresh executed odd certificate for the surplus trial, not an assumed extension of AY's even result.

## 3. Joint signed arithmetic trial

Retain AY's 640 carriers (128 low and 512 high), newly selected by leverage of proposed near-threshold directions in the odd base comparator. All carriers are actual distinct restricted circular modes. For P their synthesis, let F_P=P*J0T, G_P=P*J0R, L=F_P-G_P V, D_P=P*J0P+epsilon I. AY's variational theorem, valid for every trial Z, gives

    T*(D+B)^-1 T <= B_V - (8/35)(L*Z+Z*L-Z*D_P Z).

The initial unsplit arithmetic run completed F_P but timed out before G_P; it contributes no accepted spectral count. The replacement `az_cross_checkpoint.py` regenerates the same phase intervals, checkpoints F_P, the carrier metric, and four contiguous chunks of G_P separately, using AX's signed pairings for all 24 prime powers. `az_combine.py` checks exact index coverage, adds chunk row errors and takes the maximum column error, then forms L and the symmetric carrier metric using checked grid products. The new norm error enclosures are

    eta_F = 1323/16777216,
    eta_G = 1433/16777216,
    eta_L = 7864859317/17592186044416,
    eta_D = 1281/33554432.

The inherited AY postprocessor was executed twice on this one completed odd input. It verifies the proposed arithmetic trial, spectral Gram and residual using the exact/radix-split binary64 dot model; the numerical eigensolver is a proposal only. Both regularizations pass the strict beta<997/1000 check.

| Ridge epsilon | Trial norm bound diagnostic | Total matrix enclosure error diagnostic | Certified absorbed dimension | Odd constraints | Exact beta diagnostic |
|---|---:|---:|---:|---:|---:|
| 1/16 | 7.282517 | 0.00908306 | 4545 | 598 | 0.994481613 |
| 1/64 | 13.613829 | 0.01160977 | 4549 | 594 | 0.992502873 |

For the stronger accepted result, the exact rational beta is

    72133734561829758738096568320932443 / 72678615437072459533662337718288384 < 997/1000.

The ridge was chosen as an exploratory numerical parameter; all accepted inequalities are checked again for that trial. A global optimum over ridge, carrier choices, or subspaces is not claimed. No new carriers beyond the 640 used in the first odd trial were added in AZ.

## 4. Combined result

Odd codimension <=5143-4549=594. Combine with AY's separately certified even codimension 608 (including its one constant exclusion), since reflection parities are orthogonal for the original source form. The new conditional theorem is

    n_-(q0;D^0) <= 608+594 = 1202.

On the corresponding physical complement q0[f]>= (3/2000)||f||^2. Relative to AY's 1332, this removes 130 constraints, with the same coercivity margin. These are codimension bounds, not observed negative eigenvalue counts. AZ does not certify the remaining space positive or exhibit a negative original Weil direction.

## 5. Reproduction and provenance

The ZIP contains the AY inputs, dependencies and historical logs plus the AZ scripts and executed logs. Large derived arrays are omitted; they can be regenerated. `az_derived_sha256.txt` records the actual local checkpoints (different floating proposal implementations can generate different valid trials; reverify, rather than requiring identical floating proposal bytes across platforms). Use Python, NumPy, SciPy, Numba, under the inherited explicit binary64 integer-dot model and standard exact compiled int64 semantics within the checked range. In the extracted directory:

```
sha256sum -c array_sha256.txt
sha256sum -c aw_sha256.txt
python ay_constants.py
python ay_check.py
export OPENBLAS_NUM_THREADS=4
export NUMBA_NUM_THREADS=4
python az_prepare.py
python az_spotcheck.py
python az_solve.py
python az_base.py
python az_cross_checkpoint.py select
for x in F D G0 G1 G2 G3; do python az_cross_checkpoint.py "$x"; done
python az_combine.py
RIDGE_DENOM=16 python ay_postprocess.py 1
RIDGE_DENOM=64 python ay_postprocess.py 1
```

`az_check.py` reconciles the archived successful logs and exact rational bounds; it is not a new matrix computation. The two interrupted runs are preserved as partial logs, and `az_arithmetic.py` is the superseded unsplit arithmetic builder. The Numba compiler reports that the inner kernel is not parallelized because of its assertions; the completed integer loops and checks are unaffected. No remote CI or Lean validation.

## 6. Next gate — BA: jointly reoptimize the surplus and arithmetic trials

There is another positive correction available without introducing new carriers. In the full variational expression for r=T c, take the arithmetic trial z=B^(1/2)P Z c/21, while retaining V initially. The surplus-coordinate gradient of that expression is

    E_joint = HV-C* + R*BP Z/21 = E_V+(8/35)G_P* Z.

The quadratic cost of replacing V by V+Delta is exactly

    E_joint* Delta + Delta* E_joint + Delta* H Delta.

Choose Delta=-E_joint. Since H<=I, this change lowers the trial comparator by at least E_joint*E_joint. More generally choose Delta=-H^-1 E_joint for the exact surplus reoptimization, or use an arbitrary proposal with a residual-certified error. Thus for each previously admissible Z the new bound is

    T*(D+B)^-1 T <= B_V - (8/35)(L*Z+Z*L-Z*D_P Z) - E_joint*E_joint.

This is a valid further subtraction, even when V is only an approximate inverse trial. The gradient uses the already enclosed G_P and the surplus residual; it is not generally zero after the arithmetic trial. Alternatively, use the full exact small carrier metric P*(B+B^2)P, and its Schur complement against H, to optimize the two trial coordinates simultaneously. The displayed BA inequality is derived analytically here; its extra Gram and new spectral dimension have NOT been computed or certified. Recovering the entire remaining exceptional space is still open. Fixed-window/global Weil positivity, completed Suzuki–Connes transfer, Phase IV and RH remain OPEN.