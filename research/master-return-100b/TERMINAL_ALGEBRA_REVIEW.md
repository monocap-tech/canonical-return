# RETURN100-B — Terminal algebra review and receiving theorem assembly

Status: analytic review of the archived RETURN086 checker, including
the RETURN093-B exact-PD variant. No new numerical solve or Lean claim.
The RETURN099-B conditional physical floor is not promoted here.

## 1. Terminal trial identity: no positivity assumption on the saddle

Let K=J0+T A^-1 T* be the actual symmetric finite saddle of the
constraint/signed/restoration system (notation RETURN095-A).
Let X be the rational symmetric approximate inverse and Y the rational
trial solution used in the captured checker. Here X,Y denote the
scaled real matrices, not the integer encodings.

If either ||I-KX||_infty<1 or ||I-KX||_1<1, the Neumann series
implies K and X are invertible. In the first case use the infinity
norm; in the second use the 1-norm. Since X and the true K^-1
are symmetric, their 1-norms equal their infinity norms. Consequently
the checker's use of the smaller of the two certified residual
bounds is valid, and

    ||K^-1||_2 <= ||K^-1||_infty <= ||X||_infty/(1-rho).

The first inequality follows for a symmetric matrix from its
spectral decomposition or the general sqrt(1-norm times inf-norm)
bound. K itself need not be positive.

For U=T A^-1 m and TT=m* A^-1 m, the exact identity is

    TT-U* K^-1 U
      = TT-Y* U-U* Y+Y* K Y-E* K^-1 E,
    E=U-KY.

In particular, one cannot use the sign of the residual quadratic
form unless K is positive. Instead

    -E* K^-1 E <= ||K^-1||_2 ||E||_F² I.

The archived checker's doubled centers and absolute-radius
products bound the entries of the trial core and E at common
exact integer denominators. The symmetric trial-core error has
spectral norm at most its maximum absolute row sum. The entrywise
absolute residual enclosure bounds its Frobenius norm squared.
Thus the checker really bounds the Schur payment above, with its
positive extra term equal to the inverse cap times the squared
residual cap plus the core radius. No comparison between a
numerical eigenvalue and zero is used.

## 2. Raw H lower and exact congruence

Let Hmid, R be the symmetric raw retained center and nonnegative
radius in physical units, and w a strictly positive weight vector.
For every real y, the inequality
2|y_i y_j| <= (w_j/w_i)y_i²+(w_i/w_j)y_j² gives

    y*(H-Hmid)y >= -sum_i [(R w)_i/w_i] y_i².

This is the weighted raw diagonal majorant in the checker, with
upward integer rounding BEFORE the ill-conditioned chart
congruence. Writing C/T0 for the retained chart, the computed
Hlow is a Loewner lower bound on the actual chart energy H.
This is distinct from the scalar chart floor used for the error
norm, which RETURN099-B verifies independently.

If the final comparison
target Hlow - (4/5) PaymentUpper
has a positive definite congruence by a nonsingular V, the
same difference is positive definite. The modular full-rank
tests imply actual integer nonsingularity: a determinant
nonzero modulo the tested prime is nonzero over Q.
The alternate RETURN093-B check uses fraction-free leading
principal minors when diagonal dominance fails. Sylvester's
criterion applies since the enclosed lower matrix is symmetric.
Neither method requires invertibility of a numerically
rounded real eigenvector basis.

The factor 4/5 is NOT a chart normalization:
the physical complementary lower energy is (5/4) D in the
rephased Fourier coefficients. Its physical inverse is
(4/5) times the constrained inverse of D. The original
retained chart energy is already physical. This identifies the
normalization in the terminal payment and in all error caps.

## 3. Assembled receiving implication

For real f in the log-weighted supported form domain, let
f=p+q be the parity-preserving L2 orthogonal retained decomposition.
Write p in chart coefficients x and q in complementary coefficients a.
The mixed functional is t(x) in the same rephased coordinates.
Let D_phys=(5/4) D on the actual constraint kernel.

Assume the literal operator identification and all interval
enclosures have been validated. The actual form has

    W(f) >= H[x]+D_phys[a]+2<t(x),a>.

Let m be the full modeled source. RETURN093-B proves
||(D_phys)^(-1/2) m(x)||² <= mu H[x] with mu=49/50
even and 19/20 odd, conditional on that identification.
Here the inverse square root is understood on the kernel,
with orthogonal projection of the source onto that kernel.
RETURN085/095-A bound the three discrepancies in the same
dual norm by eps_i sqrt(H[x]). Their sum is bounded by
triangle, not by orthogonality. With r=sqrtcap+sum eps_i,

    |<t(x),a>| <= r sqrt(H[x] D_phys[a]),
    W(f) >= (1-r)(H[x]+D_phys[a]).

RETURN099-B obtains independent retained chart lower floors
and uses the crude upper mass estimate
||p||² <= (5/2)*87*||C/T0||_F² ||x||².
The bound |sinh(x/2)|<1 holds throughout I.
The complementary lower is (12/25)*(99/1000)||q||²:
the physical 5/4 factor cancels the Parseval factor.
Thus the conditional physical L2 floor is the positive
minimum of the two component floors times (1-r).
Parity decouples because reflection commutes with all
physical terms (the Mellin pole product is reflection invariant).
Orthogonality of p and q and of the two parities gives the
common conditional floor 3*10^-59 on the specified aperture.

## 4. Explicit residual obligations, not just arithmetic

The above reviews the terminal *algebra*. It does not
independently regenerate each infinite phase/auxiliary contraction
or prove the correctness of every mpmath interval primitive.
Before unconditional promotion, review at least:

* literal identification of the final piecewise diagonal, constraints,
  signed and restored rows with every assembled K/TT/U bank,
  including the boundary k=32768 and tail-only odd gauge;
* the actual scalar model error norm and its trace inputs, not merely
  the positive numeric reserves;
* all analytic premises in RETURN096-A/B and RETURN098-A/B/C,
  RETURN099-A, including the explicit-formula convention and
  source coefficient generation identity;
* permanent artifact custody, and an independent mathematical review
  if claiming a research-grade unconditional theorem.

No arbitrary-support globalization follows. A fixed-aperture
theorem cannot be rescaled to arbitrary support because the
physical prime translations and Gamma multiplier are fixed.
No native Lean theorem or RH assertion.
