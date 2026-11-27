# RETURN108-C — Global raw sampling density and zero selected endpoint

Status: internal analytic proof in the ordinary raw multiplicity-weighted
model, not a native Lean certificate or a proof of RH. Reopens
RETURN108-B at cc5b0c6c86fd7c783ef8838f0dad063853ff129f and the literal
RECON V07A, V10A, V09D sources at
206c6934536a1fcfa8ec0c661dcf93eb020e5d04. Standard unconditional zeta
facts used: discrete nontrivial zeros, finite positive analytic
multiplicities, reflection symmetry, and N(T)=O(T log T). No
explicit-formula or positivity premise is needed for the density lemma.
The fixed-aperture closed-range and continuity premises of 001G and
RETURN108-B are used only for the ensuing finite crossing conclusions.

## 1. Actual raw Hilbert space and redundant coordinates

Let Z be the set of distinct nontrivial zeta zeros, write a_rho=rho-1/2,
m_rho for the positive multiplicity, and define on complex smooth
compactly supported tests

    S(f)_rho = sqrt(m_rho) int_R exp(a_rho x) f(x) dx.

The target is l2(Z), one scalar coordinate per distinct zero.
Multiplicity weights do not introduce repeated copies. The exponents
are distinct, locally finite in C, and |Re a_rho|<=1/2.
Repeated integration by parts gives, uniformly for |Re a|<=1/2,

    |int exp(ax) f(x) dx| <= C_(f,N) (1+|Im a|)^(-N)

for all N: differentiate exp(Re(a)x)f(x) in the Fourier integral,
using compact support and uniform boundedness of the real part.
The zero-counting estimate therefore implies S(f) is square summable.
This remains true for each translate. For tau_L f(x)=f(x-L),

    S(tau_L f)_rho = exp(a_rho L) S(f)_rho.

This is precisely the V07A Mellin convention, before the unitary
Hadamard rotation of reflected pairs. In the normalized V10A
coordinate space, the critical-line second coordinates are identically
zero; remove those dummy coordinates or restrict to their orthogonal
complement K_phys. On genuine nonfixed pairs the two-dimensional
unitary rotation is onto, and on critical singletons use the one
coordinate with no duplication. Thus l2(Z) is unitarily isomorphic
to K_phys. This assertion uses the standard reflection and
multiplicity symmetry, not a claim about historical Green synthesis.

## 2. A general density lemma by a Cauchy transform

Claim: the closed linear span of S(f), for all f in C_c^infty(R;C),
is l2(Z). More strongly, translates to the right of two fixed smooth
profiles suffice (Section 3).

Suppose y in l2(Z) annihilates all samples S(f). For fixed f, write
b_rho=conj(y_rho) S(f)_rho, choosing the inner-product convention
linear in the sample. Cauchy-Schwarz gives sum |b_rho|<infty.
By translation covariance and annihilation,

    0=sum_rho b_rho exp(a_rho L), for every L>=0.

This series converges absolutely for each finite L, locally uniformly
on the real L-axis, since |Re a_rho|<=1/2.
For Re z>1/2 its Laplace transform can be interchanged with the
sum: the integral of the sum of absolute values is bounded by
sum |b_rho|/(Re z-1/2). Thus

    0= C_f(z):=sum_rho b_rho/(z-a_rho), Re z>1/2.

This Cauchy series is meromorphic on the entire complex plane.
Indeed for any compact set disjoint from the locally finite pole
set {a_rho}, the distance to that set is strictly positive,
so the sum converges normally by summability of |b_rho|.
In a disk about a single a_rho containing no other exponent,
the sum of all other terms is holomorphic by the same argument.
The residue at a_rho is exactly b_rho. The complement of a locally
finite subset of C is path-connected (a polygonal path can be
perturbed around the finitely many points in its bounded region).
The identity theorem therefore gives C_f identically zero there,
so every b_rho=0 by its residue.

For any fixed rho, one can choose a real smooth probability bump
supported sufficiently close to zero so that its transform
at a_rho differs from 1 by less than 1/2. Hence S(f)_rho!=0
for such a choice, and y_rho=0. This holds for every rho,
proving the density claim via Hilbert orthogonal complement.

No moment determinacy theorem, zero separation lower bound,
RH, uniform support sampling estimate or spectral cancellation
hypothesis was used. In particular the analytic continuation is
of the smoothed Cauchy series with l1 residues, not of an
unsmoothed synthesis series.

## 3. Two explicit profile families suffice

Use the dyadic density u_b from RETURN107-B with half-width b,
a_j=b/2^j and bilateral transform

    U_b(z)=prod_{j>=1} sinh(a_j z)/(a_j z).

The exact zero set of U_b is {2 pi i k/b : k nonzero integer}:
the first factor has precisely those zeros, subsequent factors
have subsets of those zeros, and local normal convergence plus
summable tail deviations from 1 exclude additional zeros.
At zero the value is 1. Take b=1 and b=1/sqrt(2).
Their transforms have disjoint zero sets, because their lattice
spacings have irrational ratio. Both profiles are supported
strictly inside the RETURN aperture [-5/4,5/4].

If y annihilates only the samples S(tau_L u_1) and
S(tau_L u_(1/sqrt2)) for L>=0, repeat Section 2 for each fixed
profile. Every coordinate satisfies
conj(y_rho) U_1(a_rho)=conj(y_rho) U_(1/sqrt2)(a_rho)=0
after division by sqrt(m_rho)>0. At least one transform is nonzero,
so y=0. The closed complex linear span of these two right-translated
profile families is the full raw Hilbert space. A single generic
scale avoiding the countably many critical zero ordinates also
works, but no unverified specific scale avoidance is needed.

In fact, for a target only in genuine nonfixed orbit coordinates,
one fixed dyadic profile u_1 already suffices to show that the
annihilator has zero coordinates at every noncritical zero, since
U_1 has no zeros off the imaginary axis. The second scale removes
all remaining critical-line ambiguity for full density.

## 4. Exact selected global endpoint

Let U be the unitary normalization from l2(Z) to K_phys,
and E_Pi the orthogonal projection retaining all positive
coordinates and a finite packet Pi of genuine negative channels.
The selected coefficient space K_Pi is the range of E_Pi.
Since orthogonal projection is continuous and onto its range,

    closure { E_Pi U S(f) : f in C_c^infty(R;C) } = K_Pi.

Every such test belongs to a sufficiently large supported log
carrier X_t. Its selected samples are exactly A_Pi f. Hence

    R_infty=closure(union_{t>0} A_Pi X_t)=K_Pi,
    R_infty^perp={0},
    Gamma_Pi,infty=ell P_(R_infty^perp) ell*=0.

Density is in the coefficient Hilbert norm, not in any proposed
global physical Hilbert topology. It gives no uniform norm bound
for physical preimages. If a larger redundant coefficient space
is retained, its unused critical-line dummy directions can remain
in the orthogonal complement; they have zero compression onto a
packet of genuine negative channels. Thus the selected endpoint
is still zero on the genuine packet. No statement about a
restored-background marker follows without exact identification.

By RETURN108-B, for each nonempty genuine finite packet the marker
starts at I as t decreases to zero, is norm-continuous and
strictly decreasing at positive apertures, and now has limit zero.
Every ordered eigenvalue therefore meets 1/2 at a unique finite
positive aperture and falls below it thereafter. Repeated
eigenvalues may cross simultaneously. No numerical crossing
location or uniform rate is derived. Under RH there is no such
nonempty genuine packet, so this implication is vacuous and
consistent with RH.

## 5. Explicit finite-span negative direction conditional on a packet

For a genuine negative channel i, take its unit coordinate e_i^-
in K_phys. By Section 3 density, some finite complex linear
combination f of right translates of the two dyadic profiles
satisfies

    ||U S(f)-e_i^-|| < 1/4.

Consequently ||A_+ f||<1/4 and ||A_- f||>3/4, so

    ||A_+ f||²-||A_- f||² < 1/16-9/16=-1/2.

The same estimate holds for the selected signed form retaining
a packet containing i. All selected and full samples are in l2.
After attaching the classical explicit formula with the mixed
pole normalization, this yields a compactly supported complex
physical negative direction. For the real physical form, conjugation
symmetry makes the Hermitian form real on real pairs, so
W(f,f)=W(Re f,Re f)+W(Im f,Im f); at least one real component is
negative. That last step is only asserted with the correct physical
Hermitian attachment, not for arbitrary complex signed matrices.

The argument is nonconstructive in the approximation stage:
it does not identify an off-critical zero or exhibit an actual
negative test, nor bound the needed number or size of translations.
It is a direct spectral necessity argument compatible with the
classical Weil criterion.

## 6. Disposition

RETURN108-B's global endpoint is now evaluated, not merely
rephrased as a half-bound. On every genuine selected negative
packet, the ordinary raw endpoint equals ZERO, not I/2 or
a positive operator. Thus the first-crossing exclusion program
cannot prove the half-bound by a nontrivial ordinary global range
complement. Exclusion of genuine packets is exactly the unresolved
zero-displacement question. This neither establishes the existence
of such packets nor supplies a way to prove their absence.

The two-profile density result also shows why unrestricted compact
supports destroy raw coefficient-space redundancy. It cannot be
used to assert surjectivity of any fixed supported carrier, or a
bounded global inverse. On a fixed support the range is closed
under the 001G premises and may remain proper; this is compatible
with dense increasing union and the local positivity candidate.

Recommend retiring the ordinary raw endpoint half-bound as an
independent closing supplier; it has been evaluated conditionally
on genuine packet existence. Preserve this density lemma for native
formalization and comparison with historical restored markers.
No historical source changed, no workflow activated, no independent
mathematical review or novelty claim. RH remains unresolved.
