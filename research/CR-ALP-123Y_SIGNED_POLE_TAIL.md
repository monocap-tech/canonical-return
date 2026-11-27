# CR-ALP-123Y — Signed logarithmic and oscillatory leading pole envelope

Staged module src/audits/cr_alp_123y_signed_pole_tail.py, commit 462e02786001802dc21b66a0e035dab4ad9de73c. 123W run 45929279832 executing. This derivation uses original 046 blob 554e7ab54da4545987d89aa660a4abe6c951cb63 and 123H/123M centered coefficient identities. It replaces the crude 80|fb|/v+40|T|/v+2|fb| log(1+v)/v charges, not the actual underlying comparator.

Let D0=E0+b1 S1+b0 S2/4, D=D0-sum EC_m, with the original 046 reflected coefficients. Write bt_n=ep sqrt(n)+en/sqrt(n), theta_n=pi log(n)/(2a), J0=J-1. The leading combined coefficient is

    [A+B log(j/J0)+sum_n C_n cos(j theta_n)]/v_j

where B=fb b0/sqrt(a),
A=H0+B(g0+log(J0/(4a)))-2fb D/sqrt(a)-2T sum(w_n bt_n)/sqrt(a),
C_n=-2T b0 w_n/sqrt(a).

The centered sine moments remain separate, and the H0 denominator correction adds |H0|/(4v^3). This retains cancellation among all nonoscillatory pole constants and the original signed H0.

## Gamma remainder

For z=1/4+iv/2, the difference psi(z)-log(z) is the sum-integral discrepancy of 1/(z+m). By the fundamental theorem of calculus on each unit interval its modulus is bounded by integral_0^infinity |z+t|^-2 dt <= pi/v. Hence |G-log(v/(2pi))| <= pi/v+1/(8v^2)<4/v; |R|<=1/v for v>=2048. The elementary log bound gives |G|<=v and |G-g0|<=2v there (using |g0|<6).

In 046,
  g0 c+Lambda(c-b0/(sqrt(a)v))=(b0/sqrt(a))(G+g0)/v - (b0/(2sqrt(a))) G/[v(v^2+1/4)].
Also |Es-D/v^2| is bounded by
  |D0|/(4v^4)+sum |EC_m| beta_m^2/v^4
  +|b1|(1/v+1/v^2)/v^2
  +|b0|/4*(1/v+2/v^2)/v^2.
This follows by subtracting each denominator and using the R1,S3 definitions in 046. Multiply by 2v/sqrt(a) to get the recorded inverse-square and inverse-cube gamma remainders. For EC_m beyond K=256 use beta_m^2/(beta_m^2-1/4)<=2 and exp(-2a beta_m)<=(2/5)^(4m+1); the remaining affine-geometric sum is explicit. This bound decays with v; a constant interval truncation error cannot be summed over all j.

For prime pole column subtract its exact leading -(2/sqrt(a))sum w_n(bt_n+b0 cos(j theta_n))/v. The remainder is bounded by (2/sqrt(a))sum w_n |b1|/v^2 + (1/(2sqrt(a)))sum w_n (|bt_n|+|b0|)/v^3.

## Signed L2 estimate on one parity

The leading nonoscillatory numerator squared divided by j^2 has integral from J0 to infinity equal to (A^2+2AB+2B^2)/J0. The sum over one parity is bounded by half this integral plus the total variation bound, by comparing each sample to the integral on its preceding length-2 cell. Cells start at or above J0. For a=|A|, b=|B|, the reserve is [ab+a^2+(b^2+2ab)/2+b^2/2]/J0^2, from the integral of the absolute derivative.

Partial sums of exp(i j theta) over a parity progression have modulus <=1/|sin(theta)|. Summation by parts gives an absolute bound J0^-2/|sin(theta)| when multiplied by 1/j^2 or log(j/J0)/j^2; both weights' total variation is <=J0^-2. For cosine-square diagonals use cos^2=(1+cos(2j theta))/2. For cross terms use cos u cos v=(cos(u-v)+cos(u+v))/2. The module checks a strictly positive outward lower bound for each relevant sine denominator using original 046 sincos. Divide the leading squared norm bound by (pi/(2a))^2, take upward root, and add the inverse-square/cube pole remainder norm using same parity integral bounds. Separately add all three centered sine expansion charges, NOT the old H0/pole/log charges.

## Scope

No actual trial has passed through this module yet. The module preserves signed constants and exploits oscillatory cancellation; it is not asserted to meet the required threshold. Before integrating it, run original 046 namespace selftest and verify overlap of actual leading coefficients against independent direct high-mode evaluations. No additional runner activated. 123W remains the active joint trial gate.
