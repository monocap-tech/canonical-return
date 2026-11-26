# RETURN100-C — Literal auxiliary operator identification

Review of RETURN092-B/C, RETURN088, RETURN091 and RETURN094-A.
No additional numerical solve. This note identifies the reconstructed
saddle with the specified Fourier comparator, subject to the
authenticity and validity of the scalar enclosures.

Let parity p=0,1 and let J be the complementary positive even grid.
The first 16302 coordinates are the authenticated positive even
ledger k>116 through 32766, followed by every even k>=32768.
Use A0(k)=(198000+w_k)/2000000 on the finite ledger and 87/500
on the infinite portion. In odd parity only, raise A0 to 519/1000
on 65536<k<=262144 and 1199/1000 on k>262144.
The endpoint at 32768 belongs to the constant-weight tail in the
auxiliary Gram and to the variable feature band in the mixed-source
partition; the two use the identical A=87/500. This is not an overlap
within either partition.

## Auxiliary rows and signed metric

Use rephased real grid coefficients as in RETURN096-A.
For odd positive j put h_j(k)=j/(k²-j²) in even parity and
k/(k²-j²) in odd parity. Write c=25/(16*pi²).
The odd sinh constraint uses node -c. Given the ordered constraint
nodes t_0,...,t_(n-1), the archived preconditioned constraint rows are
the rational Z chart applied to

    116^(1-p+2i) k^p / prod_(l=0)^i (k²-t_l).

The divided-difference identity in RETURN096-A proves that their
kernel is the literal retained orthogonal complement. Each such
row is bounded on l2 because its tail is O(k^(-2+p-2i)).
RETURN092-B constructs precisely these rows before computing
weighted Grams. The signed carrier rows are s_j h_j, with positive
rational scales s_j taken from RETURN063, and its signed diagonal is

    J_j = s_j² pi²/(8 v_j)       (positive weight v_j),
    J_j = -s_j² pi²/(8 d_j)      (negative weight d_j),

where weights are in physical dimensionless units, i.e. stored
millionths divided by 10^6. Thus (s_j h_j)^* J_j^-1 (s_j h_j)
is exactly the signed term 8(v_j-d_j) h_j^*h_j/pi² in the model D.
Signs of rephased odd Fourier coefficients cancel in each square.
There are 481 positive and 32 negative carriers; adding 5/6
constraint rows gives 518/519 common auxiliary rows.

For odd restoration add eight rows H=B h where B is the literal
8-by-500 rational proposal bank. RETURN094-A certifies
R=B diag(pi_upper²/(8 u_j)) B^T, with rank 8, all u_j>0 and
pi_upper>=pi. For any h the Cauchy/projection inequality gives
(Ha)^T R^-1(Ha)<=sum_j 8u_j h_j(a)²/pi². Therefore this is a
valid additional positive model energy. The restoration block
in J0 is R itself, NOT R^-1, because the model D uses R^-1.
RETURN092-C adds R to precisely the lower-right 8x8 block,
and the restoration rows are ordered after the 519 common rows.
This establishes the odd dimension 527.

## Baseline Gram and its infinite-tail evaluation

Let T denote the stacked auxiliary rows. The baseline saddle is
K0=J0+T A0^-1 T^T. RETURN092-B/C finite rows are weighted by an
outward square root of A0^-1=2000000/(198000+w_k). Their signed
base-2^18 integer limbs have per-product accumulation bound
16302*(2^18-1)^2 < 2^52. Under ordinary exact binary64 arithmetic
on these integer products and partial sums, the limb dot products
are exact. The object-integer reconstruction and absolute radius
products give an enclosure of the finite weighted Gram.

For k>=K=32768, factor every common row as k^(-2+p)
times a power series in x^-2=(K/k)^2. Each denominator contributes
a factor (1-t_l/K²*x^-2)^-1. The absolute convergence is uniform
on x>=1 because every |t_l|/K²<1. Products and rational row
combinations give the coefficients and a positive geometric
majorant on the omitted series. Each omitted tail after N terms
is bounded by E*x^(-2N). For restoration use p=1 and
sum_j B_ij=0, which cancels the constant coefficient exactly.
The series can be contracted termwise against the summable
weight k^(-4+2p); the mixed truncated/remainder and
remainder/remainder terms are bounded by the positive moment
sums at exponents increased by 2N and 4N respectively.

For power x^-a, a>1, the grid including x=1 has sum
between K/[2(a-1)]+1/2 and that value+a/(4K).
Indeed the composite trapezoidal rule is above the integral
by convexity, and its Peano error is at most (2/K)^2/8
times integral |(x^-a)''| before division by the mesh 2/K.
The latter integral is a. RETURN092-B/C use precisely this
moment formula. For the strict tail k>K the corresponding
endpoint term is -1/2; the RETURN088 cutoff convention is strict
and must not be confused with the inclusive baseline tail.

## Final reweighting

Define G_K = T P_(k>K) A0^-1 P_(k>K) T^T. These are the
strict-tail Grams independently reconstructed in RETURN088.
On the strengthened finite tail the reciprocal ratio is
(87/500)/(519/1000)=58/173; on the infinite tail it is
(87/500)/(1199/1000)=174/1199. Therefore, exactly,

    K_odd = K0 - (115/173)(G_65536-G_262144)
                - (1025/1199)G_262144.

This agrees with RETURN091 and RETURN092-C, including the
restoration block and the signed metric J0, which does not
change when A changes. All subtractions use outward signed
interval endpoints; positivity of a difference enclosure is
not presumed. The even K is the baseline common Gram plus J0.

Identical reciprocal ratios apply to the disjoint support
contributions of TT=m^T A^-1 m and U=T A^-1 m; on the odd
infinite tail the actual compensated source and its restoration
cross must be reweighted together. RETURN091 consumes the
restricted-bank signed corrections and the direct strengthened
asymptotic restoration cross from RETURN090-B. This matches
the analytic identity of RETURN098-B, rather than treating
the tail compensation as a zero functional.

## Review boundary

This identifies the auxiliary row formulas, metric, and exact
reciprocal reweighting in the reconstructed operator. It is
not an independent rerun of the scalar phase enclosures or
a new terminal matrix verification. Further review must verify
the literal TT/U finite-band phase summation and all scalar
input enclosures, including the actual pi interval in RETURN062,
and review the numerical arithmetic assumptions. In particular
the earlier terminal PASS is conditional on the mathematical
soundness of these enclosures, not a consequence of overlapping
parent matrices. RETURN099-B's conditional coercivity floor
remains unchanged. No unconditional theorem or RH.
