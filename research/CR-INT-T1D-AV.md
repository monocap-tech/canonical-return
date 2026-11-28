# CR-INT-T1D-AV — Residual-certified physical surplus inverse and low-frequency release

**Status:** executed local computational audit PASS both parities and separate low-frequency census PASS. Conditional on the original source identities AP/AM, certified inherited AQ/AS arrays, and AR's ordinary binary64 exact-integer-dot model. Not Lean-certified, independently reviewed, remote CI, or positivity of the complete fixed-window/global Weil source. No negative original direction found.

## 1. Exact variational inverse certificate

Preserve the physical interval I=(-2,2), its original two-pole-null smooth source, all 24 prime powers with their signs, the AP hybrid diagonal lower comparator (2/5 exact symbol and 3/5 AM norm), reference t=1/2 and AS's selected surplus coordinates. Let D_+=E*DE, T be the inherited high deficit synthesis, M=T*T, R_s the weighted surplus synthesis, C=T*R_s, H_s=diag(1-alpha)+R_s*R_s. AU proves B=T*D_+^-1 T=M-C H_s^-1 C*, I/24<=H_s<=I. The inherited arrays are hash checked, not rescanned.

For any rectangular trial V of size surplus by deficit and residual R=H_s V-C*, completing the square gives the exact identity

    B=M-CV-V*C*+V*H_s V-R*H_s^-1 R.

Hence B_V=M-CV-V*C*+V*H_s V is an upper comparator with 0<=B_V-B<=24||R||^2 I. This requires no trust in a proposed Cholesky factor or linear solve. Equivalently B_V=M-sym(CV)+sym(V*R), where sym(X)=(X+X*)/2. For grid approximations M0,C0,H0 satisfying norm errors eta_M,rho,eh, and any grid trial with ||V||<=v, let r0>=||H0 V-C0*||. Then ||R||<=r0+eh v+rho=:r. The symmetric grid matrix M0-sym(C0 V) is an upper comparator for B after adding (eta_M+rho v+v r+rounding)I. In addition the actual variational upper B_V differs from the exact retained B by at most 24r^2 I. The latter bound is not the full matrix enclosure used in the spectral check.

## 2. Executed full retained-surplus solve

Use all AS surplus columns. Construct the signed physical surplus-to-surplus Gram by the AQ/AR exact overlap formula, with directed rational sqrt intervals and fixed-point row enclosures; add diag(1-alpha). Check symmetry exactly. S=2^26 grid for M,C,H, U=2^20 for V, 2^24 eigenproposal grid. Dense Cholesky solve proposes V only. The full residual H0 V-C0* and cross product C0 V are checked in blocks with integer operands represented as ordinary IEEE binary64; Cauchy bounds on sums of absolute products <2^53 ensure exact integer dots, and grid conversion and int64 overflow bounds are checked. The final symmetric cross-product floor adds at most m/S. The signed spectral subspace is proposed by eigensolve, then independently checked with AR's exact Gram and residual inequalities. All acceptance comparisons are exact fractions. No solver accuracy is assumed.

| parity | surplus dimension | deficit dimension | H norm error (diagnostic) | residual upper diagnostic | absorbed dimension | constraints incl. 512 | beta diagnostic |
|---|---:|---:|---:|---:|---:|---:|---:|
| even | 12560 | 4628 | .0000949055 | .0023470692 | 4132 | 1008 | .9958803411 |
| odd | 12751 | 4691 | .0001900197 | .0035138548 | 4185 | 1018 | .9918461431 |

Both final audit processes exit 0. The exact log checker proves beta<997/1000, residual r<2348/10^6 even and <3514/10^6 odd, and 24r^2<1/7500 even and <3/10000 odd. Conditional inverse-only negative index <=2026 (AU 2245). All diagnostic decimals are informational; exact comparisons use the logged fractions.

## 3. Independently release 119 blanket low exclusions

AQ through AU excluded every circle coordinate 0<=j<1024, 512 in each parity (frequency j*pi/8). AV checks these coordinates individually; this does not require their physical overlap with any high modes. Use the exact digamma expansion with x=j*pi/16 and lower pi and m_Gamma(0) from AM. Sum the first 32 positive terms x^2/[(k+1/4)((k+1/4)^2+x^2)] exactly. The remaining decreasing tail is bounded below by the integral (1/2)log(1+(x/(129/4))^2). At x=0 this tail is zero. Round the log argument downward on a 10^9 grid and evaluate by the rigorous AM logarithm series. All 24 signed prime phases are bounded by AQ's certified Taylor and recurrence arithmetic (direct low-index recurrence, with the inherited conservative 1/150 arithmetic symbol slack). For each j certify the hybrid diagonal d_j >= t or leave it excluded. No use of the high-index asymptotic floor at these low frequencies.

The executed low census releases 59 even and 60 odd modes, leaving 453 and 452 low constraints. Its minimum released surplus lower bound is strictly positive (diagnostic .0005393562). The released coordinates have d_j>=t; their surplus above t can simply be discarded without weakening AS's existing lower comparator on the enlarged admissible subspace. The high deficit and surplus arrays, high inverse, and accepted high-coordinate subspaces are unchanged. No low negative deficit is absorbed, and no prior Ritz dimensions are subtracted.

Combined conditional codimension <= (1008-59)+(1018-60)=949+958=1907. For f on that certified complement, T*f lies in the certified subspace, so ||T*f||^2=<f,T(T*f)> <= <D_+f,f>^(1/2)<B T*f,T*f>^(1/2) gives ||T*f||^2<=beta<D_+f,f>. Thus q0[f]>=3/2000 <D_+f,f> >=3/2000 ||f||^2. Conditional negative-index upper bound n_-(q0;D^0)<=1907. This counts possible negative directions, not observed negative eigenvalues.

## 4. Reproduction

In the bundle directory, Python 3/NumPy/SciPy:

    sha256sum -c array_sha256.txt
    OPENBLAS_NUM_THREADS=4 python av_inverse.py 0 > av_even.log; echo $? > av_even.exit
    OPENBLAS_NUM_THREADS=4 python av_inverse.py 1 > av_odd.log; echo $? > av_odd.exit
    OPENBLAS_NUM_THREADS=1 python av_low.py > av_low.log; echo $? > av_low.exit
    python av_check.py

The full solve needs several GB of RAM and a temporary grid file (removed after residual check). The local runs above actually completed. The log checker checks these particular logs and is not a substitute for rerunning the auditors. The inherited arrays can be regenerated by aq_scan.py and as_surplus.py. No kernel or remote CI claim.

## 5. Next gate — AW: remove the remaining low-coordinate blanket constraints with the original source energy

905 of the remaining 1907 constraints are low circle coordinates whose hybrid lower diagonal is not certified above t. The retained high deficit/surplus inverse is now solved to a small certified variational residual. Continuing exterior iterations cannot recover the discarded positive source terms. AW should retain the low deficits explicitly, including the normalized constant mode and exact two pole moments, or construct a joint low block using the original boundary-restored Weil matrix and its mixed terms. A full exact inverse of the retained hybrid comparator is still only a lower comparison: it discards positive Gamma image and arithmetic slack, and a negative comparator direction would not falsify original Weil positivity. Full fixed-window/global positivity, completed Suzuki–Connes transfer, Phase IV and RH remain OPEN.