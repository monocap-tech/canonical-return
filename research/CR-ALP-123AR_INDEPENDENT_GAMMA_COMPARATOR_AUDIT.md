# CR-ALP-123AR — Independent full-frequency Gamma comparator audit

Research-only; no fixed-collar or RH promotion. Source 068 band payload blob 1a7fe8629daa8b8e6cec185a79b6746d6b28975d. Original resolvent derivation 122G and 066/068 centered-boundary estimate reviewed. No change of actual comparator or trial coordinates.

## Direct symbol lower

For z=1/4+iv/2, v>0, the digamma series gives psi(z)=lim_N(log N-sum_{m=0}^{N-1}1/(z+m)). The integral over [0,N] of 1/(z+t) is log(z+N)-log z with the principal branch on Re z>0. Subtracting each left rectangle from its integral and applying the fundamental theorem of calculus gives modulus at most integral_0^infinity |z+t|^-2 dt <= pi/v. Thus Re psi(z)>=log |z|-pi/v>=log(v/2)-pi/v. Machin and positive exponential Taylor sums certify pi<22/7 and exp(46/25)>44/7>2pi. Hence g(v)>=log v-46/25-22/(7v). The right side increases for v>0; for v>=k>=120 it exceeds log k-37/20-2/k because 1/100-8/(7k)>0. No historical staircase bound needed.

## All-resolvent comparison

For each parity, the exact pole-null identity annihilates f_(1/2)(v)=v/(v^2+1/4) on high sine coefficients. Subtract alpha_beta f_(1/2) from f_beta with alpha_beta=(120^2+1/4)/(120^2+beta^2). The resulting factor is f_beta times (beta^2-1/4)/(120^2+beta^2) times (v^2-120^2)/(v^2+1/4). Apply weighted Cauchy-Schwarz with positive weights f_beta(v)/sqrt(v), keeping the squared centering factors on the second factor. For beta<32 use the independently recomputed small reserve 6631424/5028336225; for beta>=32.5 the unimodal lattice bound follows by integrating each monotone side and allowing one peak rectangle. The integral of sqrt(v)/(beta^2+v^2) on the positive line is pi/sqrt(2 beta). The remaining beta summand is decreasing, and its full step-two integral yields I=pi/(2 sqrt(2)); the tail beyond split B is bounded by v^(3/2)/(3(B-2)^(3/2)). Thus each saved split is a valid bound using the original 068 inequality; replacing the square roots by new directed rational endpoints at 10^18 yields no larger penalty than the saved one. The last band uses the full universal penalty and no frequency cutoff.

The exact BigInt execution recomputed 118 contiguous bands, checked every saved lower endpoint against the direct symbol minus independently evaluated penalty, including the infinite final band. 32-term positive atanh logs; Machin alternating bounds and positive Taylor exponentials. No floating acceptance (integer bit-length initializes square-root iteration only).

| Quantity | downward rounded to 1e-9 |
|---|---:|
| Minimum saved gamma lower | 2.890358630 |
| Minimum independently derived gamma lower | 2.896809385 |
| Minimum direct-symbol improvement over inherited lower expression | 0.000632319 |
| Minimum saved gamma from 2049 onward | 4.447938456 |

The capped comparator remains unchanged. With the historical prime operator norm <=21/8, its uniform floor is at least 51/200 and its post-2048 diagonal-minus-prime floor is exactly 1.822938456. Prime norm is an inherited premise here; original coefficient/domain and tail-envelope proof acceptance remain distinct. No claim of native Lean certification or global positivity.

## Custody

Source src/audits/cr_alp_123ar_independent_gamma_bands.js, commit 9513dcb42cd695be78323dce67775b1f8fe405e3, blob 3eaf29b481dda966a9c93b62f1b51d29eff3fbd8. Receipt custody/cr-alp-123ar-gamma-bands/result.json, commit be2c2f6621c71fcc1b64450648a2378a9ee23b90, blob 603191708a0c265b3508112a1579f88559c50ef3. Both exact remote content readbacks matched. Replay with node src/audits/cr_alp_123ar_independent_gamma_bands.js src/audits/CR-ALP-068_pole_centered_gamma_result.json /tmp/ar-result.json. Execution used the identical audit function in a JS BigInt runtime; standalone node launcher has not been executed here.

Separate AQ execution failed after cosh coefficient comparisons due to mat name shadowing. No dependence of this audit on AQ matrix or response outputs.
