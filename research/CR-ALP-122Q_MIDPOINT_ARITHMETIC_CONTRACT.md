# CR-ALP-122Q — Physical midpoint arithmetic contract readback

Research-only. The quantitative receiver 45917791950 passed, durable 6e04c83, conditional physical constants 1e-38 odd and 1e-34 even. No new runner or positivity promotion. Reopened actual frozen 058 validated.cpp db1c9bac0c05cd82065d42fb8b5386245cda8d5b, 058 VALIDATION.md f45e7cc9c523a41525842728626be7fcd896fb14, and 088 response.cpp 7380bb7084855f0b043f202aa412d2fbfb859794. Distinguish the ideal quadrature, interval arithmetic, and midpoint fast-path obligations. The first two have exact mesh and bracket evidence, and the third requires the explicit classical forward-error argument under the hardware/library assumptions below.

## Arithmetic boundary

Assume round-to-nearest IEEE binary128 basic operations without fast-math/reassociation (including subnormal support), exact nextafterq adjacency, and correctly rounded 192-bit MPFR with the declared ABI and rounding codes, with directed conversion to binary128. These are trusted implementation premises rather than Lean theorems. For normal results a conservative relative unit u=2^-112 suffices; for subnormal results use absolute 2^-16494. Directed interval basic operations enclose by one nextafter on each side; the midpoint trigonometric enclosure uses the radius and the 1-Lipschitz property. The decimal grid multiplication is outward before floor/ceil, provided the result fits int64. The native response code bounds its printed G/E/Y2/Yb/y/pole/X2 before printing; the original moment code also requires a range argument for its corr, M and penalty outputs (see below). The runtime checks in 122K are not by themselves proof of forward error.

## Midpoint budget check

Let exact arguments and coefficients be chosen within the validated input intervals (after first projection, likewise for the updated interval pieces). The numerical midpoint error is additional to the explicitly accumulated coefficient radii. In eval the restarts every 16 sine indices use angle error <=850*(2.2e-32+(10/7)*1e-31)+1e-29<2e-28. The step angle error is <4e-31. Correctly rounded MPFR followed by binary128 rounding incurs at most two binary128 ulps. Over at most 15 rotations a Euclidean vector error <1e-27 is sufficient (exact rotation norm one; per-step rounding <20u). Times the checked sum of sine/cosine piece magnitudes <1e17 this contributes <1e-10. The normalized basis padding 1e-26 covers division by sqrt(a) and the same phase budget.

For the 128 low exponentials, the exact argument perturbation is controlled by ti.rad<1e-31, d.rad<1e-32 and the binary128 literal 2.2. All recurrence factors have magnitude at most one on the physical interior. A per-exponential absolute budget 1e-27 covers the initial errors, at most 127 multiplications and phase displacement; the reflected sum times the verified low amplitude sum <1e14 costs <2e-13. For grouped Poisson terms the real argument error is <4e-20 by c.rad<1e-20, c<2^34 and the time radius. The real Poisson derivative is bounded by 2; reflected amplitude sum <1e3 gives <2e-16. For x<=1000 the recurrence has at most 39 multiplications/divisions, each relative error at most u for normal results; allowing 200u times term magnitude <=1 is ample. In the skipped x>1000 branch, even exact x>=999 gives e^-x x^j/j! <= (10/27)^999 1000^39/39! <1e-300 for 0<=j<=39: monotone in x for x>=999 and increasing in j there. e>27/10 follows from its first six Taylor terms. Subnormal intermediate error amplification <=1000^40*2^-16494<1e-300. The bounded resonant and ordinary exponential pieces have argument error <1e-30, magnitudes <4 for individual exp(+-t/2), and piece coefficient sum <=2e17. Allow <1e-11 for these terms including resonant arithmetic. The remaining fewer than 5000 scalar operations on partial magnitudes <2e18 cost <2e-12 in the stated classical budget. Even allocating 1e-9 for the trigonometric term and rounding the remaining allowances upward gives a total <1e-8. The checked coefficient-radius reserve in midpointize is separate. This readback follows the original VALIDATION.md estimates, not a machine-checked operation-by-operation rounding proof.

Exact rational side checks executed in JS BigInt:
[
  [
    "rotation arithmetic reserve",
    true
  ],
  [
    "skipped Poisson tail at x>=999",
    true
  ],
  [
    "binary128 subnormal amplification",
    true
  ],
  [
    "output grid range",
    true
  ],
  [
    "point error sum",
    true
  ]
]

## Output range for original moments

The original first-moment code has a runtime bound G.h<1 after quadrature padding, before it evaluates shifted correlations. Conditional on the preceding enclosure, ||rho||^2<1. On each reflected half-overlap the positive ideal quadrature sum of either squared proxy factor is <=2||rho||^2+1e-12, by the same analytic-piece Cauchy estimate for its square, not a pointwise global amplitude bound. Interval eval endpoints differ from the exact value by at most 0.02 (plus negligible outward rounding), since midpointize checks radius <0.01. Discrete weighted Cauchy-Schwarz, total reflected overlap length <=2.2, gives absolute correlation endpoint bound <3 including the stated 1e-12 padding and ordinary interval rounding. Hence |M|<24 with paired weight sum <8 and |penalty|<(2.88+24)/(0.255*5.505)<20, far below int64 grid capacity. This is a classical conditional range argument; an explicit pre-conversion assertion in out would be a stronger executable guard. The preserved output logs were already replayed, but successful casts alone do not prove range.

## Status

The prior 058 classical forward-error premise is now explicitly separated from the exact arithmetic receiver. Still not a formal IEEE/MPFR verification, exhaustive independent reintegration or proof that all continuum coefficient identities and generators attach to the exact same source in Lean. The very small physical coercivity bound reflects the large conditioned-coordinate norm, not a numerical instability of the final exact congruence. Do not infer unrestricted positivity or RH.
