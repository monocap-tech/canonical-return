# RETURN108-G — Reflected duplicate quotient and marker endpoint audit

Status: exact internal algebra for the canonical duplicate convention;
historical executable normalization still not recovered. No native Lean
verification, no RH promotion. Refines 108-C/F: zero endpoint applies
to genuine active coordinates, not arbitrary redundant compressions.

## 1. Source audit

Reopened monocap-tech/weil main:
- WeilDefect/PairGeometry.lean blob 5361e8384880fcd799fbd62151b65be7a11f3cf8.
  pairEigenEquiv is an unnormalized linear equivalence, inverse divides
  sums/differences by 2; it is NOT a unitary normalization theorem.
- docs/QUARTET_CHANNEL_RESIDUE_STRUCTURE.md blob
  6a6dad007688e8415f31aa660bb19c5a6644dff5, sections 1, 2, 4:
  canonical basis (e_gamma +/- e_conj(gamma))/sqrt(2), one positive
  fixed-point coordinate, two genuine negatives per simple quartet,
  same-frequency multiplicity null quotient.
- docs/IMPORTED_SOURCE_PINS.md: finite inertia/multiplicity external
  premises explicitly distinguished from internal native Green metric.
- Library "Pasted text.txt", archived actual-zero mission, statements
  actual_zero_negative_form_obstruction and negative_margin_forbids_marker_lower:
  these expressly retain reflection duplication, original multiplicities,
  pair normalization and actual selected packet. No source body for the
  signed actors is present in the retrieved snippets.
- Library "Pasted text(4).txt" raw source pairing as in 108-E.
The archived Git pin f02a526d7c7fb8a43356376a2389d6f8afc08e2a returned
422 in the previous weil tree lookup and again in the weil-lab commit
lookup. Search for literal signedColumn/positiveActor/negativeActor
did not recover executable bodies. Do not promote a historical source
identity from the mathematical reconstruction below.

## 2. Canonical active and duplicate coordinate realizations

Take a complex Hilbert physical space and any raw columns c_rho,
c_Rrho on a nonfixed reflected orbit, with m_rho already included
through sqrt(m_rho). Let c1,c2 abbreviate the two. The canonical
orbit columns are p=(c1+c2)/sqrt(2), n=(c1-c2)/sqrt(2).
Their signed covariance is p p* - n n* = c1 c2* + c2 c1*,
where a b* denotes the rank-one map x -> a <b,x>.
Their unsigned covariance is c1 c1* + c2 c2*.

A full reflected-index duplicate convention that preserves both
covariances uses positive columns p/sqrt(2),p/sqrt(2) and negative
columns n/sqrt(2),-n/sqrt(2), i.e. the raw sum/difference divided
by 2, not sqrt(2). At a fixed reflection point use the original
column c and zero negative column (or omit that coordinate).
These definitions give the exact original per-orbit signed and
unsigned covariance. Multiplicity is already included in c:
do not multiply a second time. If m repeated identical unweighted
columns are retained, the active isometric embedding of one
weighted coordinate is the constant vector 1/sqrt(m) per copy,
whose synthesis is sqrt(m)c_unweighted; orthogonal zero-sum
vectors are null. Equal multiplicities under reflection ensure
the two reductions commute.

Within the specific symmetric sum/difference scalar convention,
the normalization is forced by equality of the orbit covariance:
one pair representative needs squared coefficient 1/2 on the
raw sum/difference, two copies need squared coefficient 1/4.
This is not a uniqueness theorem among arbitrary indefinite
changes of basis. It presupposes the stated canonical channel
geometry and positive metric.

## 3. Exact coisometric quotient

Let K_+ and K_- be genuine orbit-reduced coefficient spaces.
Let K_+^d,K_-^d be duplicated spaces. Define isometries V_+,V_-
from active into duplicated coordinates, per nonfixed orbit

    V_+ u = (u,u)/sqrt(2),
    V_- u = (u,-u)/sqrt(2).

On fixed points V_+ is the identity, with no active negative
coordinate; any dummy negative coordinate is in the complement.
These direct-sum maps are isometries. If P,M are canonical positive
and selected negative syntheses then

    P_d=P V_+*,  M_d=M V_-*

on reflection-complete selected packets. This follows columnwise
and extends by boundedness; hence P_d P_d*=P P* and M_d M_d*=M M*.
For arbitrary reflection-complete complementary tails the same
calculation preserves the covariance and relative tail bound.
All physical signs and screening comparisons transport.

For a reflection-complete selected duplicate space let V denote
its negative embedding. Given any uniformly positive covariance
Q (ordinary or restored), set C=(I+M* Q^-1 M)^-1. The exact
duplicate marker is

    C_d=(I+V M* Q^-1 M V*)^-1
       =(I-V V*) + V C V*.

Proof: split into Ran V and ker V*; on the first the operator is
unitarily conjugate to I+M*Q^-1 M, on the second it is identity.
Thus the active marker is V* C_d V=C. Its lower half-bound is
equivalent to C_d>=I/2, since the redundant block is I.
The ordinary gap endpoint has the same decomposition by taking
the fixed-window regularization limit. Under 108-C/F, conditional
on nonempty genuine negative packet and correct source attachment,
the global active limit is zero and the full duplicate marker
limit is I-VV*, NOT zero in norm when redundancy is present.
Its active eigenvalues go to zero, its redundant ones equal 1
at every stage. This also applies to the simultaneous restored
limit under the relative-tail admissibility premises.

## 4. Incomplete duplicate packets: compression versus re-selection

On one reflected negative pair the limiting redundant block is

    I-VV* = (1/2) [[1,1],[1,1]].

Compressing this FULL two-coordinate marker onto one of its
duplicate coordinates gives precisely 1/2. At finite parameters
the one-coordinate compression is (1+C)/2 for a scalar active
marker C. This stays above 1/2 whenever C>0, despite C eventually
crossing its own half-threshold. In a larger packet the diagonal
entry is (1+<C e_i,e_i>)/2 with same consequence.
This is an exact warning, not a claim that the original
implementation used an incomplete packet.

Re-selecting a single duplicate actor BEFORE forming the
Picard inverse is different. For one orbit with canonical
cost k=M*Q^-1 M>=0 scalar, the new selected synthesis is
M/sqrt(2), so its marker is (1+k/2)^-1=2C/(1+C), where C=1/(1+k).
The full duplicate marker compression is (1+C)/2 instead;
their difference is (1-C)^2/[2(1+C)]>=0. The re-selected
marker tends to zero as k->infty; it is not the half endpoint
of the compression. Neither operation can be substituted for
the other by a change of notation. For finite selected packets
use reflection-complete custody or the active quotient before
interpreting the original half-marker criterion.

## 5. Historical adjudication and boundary

The recovered historical public source establishes canonical
mathematical orbit geometry; the earlier mission statement records
full critical-strip indexing and preservation of duplication and
normalization. The precise original signed actor code and whether
every historical selected packet was reflection complete remain
unverified. The theorem above specifies the exact audit test:
compare actual columns with the two conventions, check V*V=I,
P_d=P V_+*, M_d=M V_-*, and verify selected/complement packet
custody before using any marker endpoint. The unnormalized
pairEigenEquiv alone does not provide these identities.

Important: the spurious 1/2 is a possible redundant-coordinate
compression artifact, NOT an actual historical error established
by the available source. A reflection-complete full duplicate
marker still has a genuine negative direction below half
conditional on an off-critical channel and the global limit.
The active ordinary/restored endpoint obstruction survives
multiplicity and reflection duplication. No claim that actual
off-critical zeros exist. No new runner. Next audit should
locate the original signed actors and packet predicate or
record source custody as blocked and preserve the conditional
congruence rather than equating these models without evidence.
RH remains unresolved.
