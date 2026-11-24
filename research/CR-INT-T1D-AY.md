# CR-INT-T1D-AY — Joint deficit inverse and signed arithmetic-metric trials

**Status:** new analytical comparison, conditional on the original source identities and the AM/AQ/AS/AV/AW/AX inputs. Fixed-window smooth double Mellin-pole-null Weil form on I=(-2,2). All 24 original prime powers retain their original signs. No Lean, remote CI, independent mathematical review, or proof of the full fixed-window or global positivity statement.

## 1. Joint source comparator

Use AX's reference energy t=1/2, lambda=399/50, normalized arithmetic slack B=(6/5)(lambda I−A), and J0=(lambda I−A)/4. AM gives ||A||<lambda and therefore 0<=J0<=lambda/2 I<4I and I+B<=21I. B=(24/5)J0. The positive circular Gamma image correction remains unconsumed. Let T synthesize the retained normalized deficits across **both** low and high coordinates of one reflection parity, and let R synthesize both sets of positive surplus coordinates with sqrt(alpha) weights. Then D0=diag(1−alpha)>0, D=I+R D0^-1 R*, M=T*T, C=T*R, H=D0+R*R. On the physical source with the even circular constant excluded,

    q0[f] >= t(<(D+B)f,f> − ||T*f||²)+G_sigma[f].

The exact input census uses the inherited SHA-256-verified AQ/AS arrays and AW low-frequency floor. There are 5080=452+4628 even and 5143=452+4691 odd deficit coordinates, with 12619=59+12560 and 12811=60+12751 surplus coordinates. The constant is the only remaining blanket exclusion, in the even sector. All positive low-frequency entries are retained as surplus. No new all-frequency scan was performed.

## 2. Joint variational trials

For any physical r and arbitrary surplus and arithmetic trial coordinates v,z, the inverse variational principle yields

    <r,(D+B)^-1 r> <= ||r−Rv−B^(1/2)z||² + <D0 v,v> + ||z||².

This follows by minimizing the right hand side over both coordinates; no inversion of B is needed. Take r=T c and v=V c for any matrix V. Let W=T−RV and

    B_V = M−CV−V*C*+V*HV.

For any finite synthesis P of distinct normalized parity-circle modes restricted to I, Parseval gives PP*<=I. Define U=J0 P. Then UU*<=J0²<=4J0<=B and U*U<=16I. Applying the variational principle to the weaker comparator D+UU* and choosing the arithmetic coordinate U*W c/17 gives

    T*(D+B)^-1 T <= B_V − (U*W)*(U*W)/17.       (AY.1)

This is a joint comparator: its additional term has low, high, and mixed blocks. No separate high-side certificate or mixed norm factor of 16 is needed when the joint spectral sign is certified.

## 3. Signed arithmetic-metric trial

A more flexible trial uses B directly. Let L=P*J0 W and D_P=P*J0P+epsilon I, epsilon>0. Take z=B^(1/2) P Z c /21 for an arbitrary finite trial matrix Z. Since B²+B<=21B, expanding the variational expression shows

    T*(D+B)^-1 T <= B_V − (8/35)(L*Z+Z*L−Z*D_P Z).    (AY.2)

Indeed, the physical arithmetic residual is W−BPZ/21. The cross gain is 2 Re <P*BW c,Zc>/21, and the quadratic cost is at most <P*BP Zc,Zc>/21. Adding epsilon I only weakens the gain. The numerical maximizer is Z=D_P^-1 L, which would give subtraction (8/35)L*D_P^-1 L. No trusted inverse or lower eigenvalue of the unregularized carrier metric is necessary: any rational-grid Z has the stated variational validity. The scalar and metric trials are separate sufficient upper comparators; a uniform strict Loewner ordering between their corrections is not asserted.

The carriers need not be the original exceptional coefficient basis. Their numerical selection may target eigenvectors of an approximate base comparator near the positivity threshold, as long as the final signed cross matrices and trial products are rigorously enclosed. Here P consists of selected distinct actual restricted circle modes, and its arithmetic pairings use all 24 original prime powers via AX's exact shift formula.

## 4. Numerical verification design

`ay_constants.py` verifies the scalar constants and the complete joint coordinate census exactly. `ay_joint.py` constructs a directed-grid approximation of M,C,H, chooses a Cholesky proposal for V, verifies its integer-grid residual under AR's explicit ordinary binary64 exact integer-dot model, and bounds B_V. It then applies the scalar trial (AY.1). `ay_metric.py` uses the same joint base calculation but instead selects 512 high and 128 low carriers by leverage of proposed base near-threshold eigenvectors, encloses F=P*J0T, G=P*J0R, and the regularized carrier metric, and forms L=F−GV. It chooses a proposed Z and verifies (AY.2) by directed product error propagation. Radix splitting of integer products is used if a direct binary64 dot exceeds the proved exact range; each radix product is independently bounded and final recombination is checked to fit int64. The final joint spectral test is AR's exact Gram and residual certificate on a proposed eigensubspace. No acceptance from unverified proposals.

For the base trial, with R_V=HV−C*, B_V=M−sym(CV)+sym(V*R_V). Thus the audit encloses the base trial using the cross-matrix rounding error and a certified residual norm. The sign test concerns the upper comparator (not a computed spectrum of q0), and must pass separately from the residual and arithmetic-factor checks.

If a d-dimensional subspace of the m-dimensional deficit coefficient space has upper spectral bound beta<1 for the joint comparator, impose T*f in that subspace. Cauchy–Schwarz in the D+B energy metric gives ||T*f||² <= beta <(D+B)f,f>. Hence q0[f] >= t(1−beta)||f||². Codimension at most m−d, plus one for the even constant. Parity results can be combined with an independently inherited certificate in the other parity, but an unfinished run contributes no new dimension reduction.

## 5. Executed even certificate and inherited odd frontier

The initial scalar trial `ay_joint.py 0` completed the joint base inverse and 384 arithmetic cross factors but rejected a final product at the binary64 integer-dot Cauchy range check. The first signed metric version `ay_metric.py 0` completed the joint matrix with 640 selected carriers but used a loose row/column trial norm (~160.78) and subsequently rejected its final residual product. Neither failed run supports an accepted spectral count. The final `ay_refined.py 0` uses a ridge of 1/16, an exact integer Gram bound for the proposed arithmetic-trial norm, and radix splitting for both large trial products and spectral residual checks. It completed with exit code 0.

| Item | Even result |
|---|---:|
| Deficit coordinates | 5080 |
| Surplus coordinates | 12619 |
| Arithmetic carriers (high + low) | 512 + 128 |
| Joint base surplus residual | 243484475581/140737488355328 < 1731/10^6 |
| Base surplus trial gap to exact retained-surplus inverse | <1/13000 |
| Refined arithmetic trial norm upper diagnostic | 7.489413261413574 |
| Total matrix enclosure error diagnostic | 0.019966297523889797 |
| Absorbed joint dimension | 4473 |
| Exact spectral bound | 4305399352917154040261679110710943/4326516849030769960482617363005440 < 997/1000 |
| Conditional exceptional constraints | 5080−4473+1=608 |

All quoted decimals are diagnostics, not exact endpoints. The arithmetic matrix uses all 24 prime-power signed pairings; the trial is restricted to 640 carriers, so the exact full arithmetic inverse is not evaluated. The base residual bound implies the gap <=24 times its square, as in AV; this is for the retained surplus inverse alone, not the full arithmetic inverse. The accepted spectral test verifies the joint trial upper comparator, including all its mixed blocks.

The odd AY run has **not been executed**. Retain AX's separately accepted odd 724-constraint certificate, whose exact log comparisons were rechecked by `ay_check.py`; the original odd matrix certificate was not regenerated. The two parity forms are orthogonal, and each supplies the same physical coercivity margin 3/2000 on its respective complement. Thus the conditional negative index is at most 608+724=1332, a reduction of 111 from AX's 1443. No claim that 1332 negative directions exist or that these exhaust actual negative directions. No full fixed-window positivity or negative original witness.

## 6. Reproduction and next gate — AZ

Use Python/NumPy/SciPy, under the explicitly stated AR ordinary binary64 exact integer-dot arithmetic assumptions. In the extracted directory verify `sha256sum -c array_sha256.txt` and `sha256sum -c aw_sha256.txt`; run `python ay_constants.py`, `python av_check.py`, `python ay_check.py`. To regenerate the accepted even computation run `OPENBLAS_NUM_THREADS=3 python ay_refined.py 0` (the script writes large intermediate grid checkpoints). The included failed logs are provenance, not evidence of a successful certificate. `ay_joint.py` and `ay_metric.py` are exploratory versions; use the refined script for acceptance. To reuse a builder-generated checkpoint `ay_metric_input_0.npz` for alternate arithmetic regularization, `ay_postprocess.py` implements the trial and spectral checks; that separate postprocessor was written but not executed. Large derived matrices are omitted from the reproducibility archive and can be regenerated from the included scripts. Inherited AX odd log checks are not a substitute for independent reproduction of AX's full arithmetic enclosure (see its separate bundle).

**Next AZ:** run the joint refined arithmetic-metric certificate on the odd parity, then choose additional high-side carriers by actual residual/eigen-subspace information. Reuse checkpoints for a further jointly optimized trial, rather than restarting the large retained-surplus solve for every parameter. A full arithmetic-inverse or source-equivalent coercivity theorem must still address the remaining exceptional space. Fixed-window/global positivity, completed Suzuki–Connes transport, Phase IV and RH remain OPEN.