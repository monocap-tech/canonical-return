# RETURN108-H — Surviving full-divisor normalization and arbitrary packet weights

Status: actual source recovery in a surviving Weil-Lab lineage, plus internal
analytic deductions. Not a recovery of the missing ConnesGreen historical
implementation; no new native build or strict named axiom audit.

## 1. Search custody and positive discovery

Organization code search for actual_zero_half_bound_iff_finite_restoration,
actualGreenSource, CanonicalGreenNeutralShell, ConnesGreen and
actual_zero_diagonal_marker_recovery returned no results. Commit searches
Green/RG0 in weil and Connes in weil-lab returned no results.
Branch search for connes/rg0 in weil, weil-lab, canonical-return was empty.
Branch search for green found IV091 and IV092 in weil-lab. Reopened IV091
tree and pinned ref f2e8ec89ae4bf8687a652b299b730d6d1baa43f7.
The old f02a526d7c7fb8a43356376a2389d6f8afc08e2a remains unavailable.
This is not an exhaustive audit of every historical ref or archive.

All following paths under WeilDefect/Arithmetic/ at the IV091 pin:
- ActualZetaZeroCoordinates.lean b902f7916d1e1fe76562e2ce9d0314b3522dc9d2:
  NeutralActualZetaZeroPoint = {rho | zeta(rho)=0 and 0<Re rho<1},
  ordinate=-i(rho-1/2).
- ActualZetaMultiplicity.lean 58b65e72c5d215171326eb465343e3cae9833ad2:
  analyticOrderNatAt, positivity, Sigma rho (Fin multiplicity).
- ActualZetaMultiplicitySymmetry.lean 85971f8cb856232225a804a969e9f7f42239bc1f:
  reflection and conjugation multiplicity equality, divisor involution
  preserving copy label and conjugating ordinate.
- ActualZetaSourceGraph.lean 5a19c9649315241b9e8cdf0dd3f10821cdd8db27:
  actual closed graph of physical log vector and both l2 sample families,
  coordinate equalities, unnormalized energy factor two.
- ActualZetaSelectedBackground.lean 76e24c6810241c153134c525f5ab75177015a53b:
  arbitrary finite coordinate selection (not an orbit-completeness
  predicate), exact complement, second normalization 1/sqrt2,
  selected normalized negative energy = (1/2) sum selected raw P/N
  negative source energies.
- ActualZetaRawSamplesToFullSourceGraph.lean
  66ac4401448f08af37cde97de012609bd190bec9:
  both raw sample square summabilities imply exact same-vector lift;
  the general lift is NOT limited to the separate real-odd corollary.
- ActualZetaFullSourceLiftOfSquareSummable.lean
  7a30f2c73992040de0dbf3d184f63334b6efcb2f:
  explicit full-source construction and physical identity.
- ActualZetaFullGraphRawDivisorQuadratic.lean
  5070335202dfdb22de66f01b254571c80e24c7dd:
  original normalized P/N quadratic = real raw divisor cross term;
  literature RHS on actual graph vector with C2 correlation.
- ActualZetaGreenEnergy.lean 874c9b8b103fddaded5887935d910d6e4f8d06:
  actual Green value-plus-gradient energy and summability.
- ActualZetaGreenPacketDensity.lean 2325cfc5812ee50a993de4e4a382b7cc0e3d532f:
  closed span equals Green graph closure, NOT whole graph.
- ActualZetaGreenCopyObservability.lean 8cc0d76d5f9d8bd4d763fa00c6db6e61dec14c37:
  equal ordinate copies give equal graph columns and observations.

These are executable Lean source definitions and proof bodies inspected,
not a fresh build or verification receipt. IV091 branch is used as an
immutable source witness, not current IV106 certification status.

## 2. Exact normalization recovered for this source graph

Write z=-i(rho-1/2), f the common physical vector,
v=windowEvaluation(a,z,f), u=windowEvaluation(a,conj(z),f).
The source code's first P/N mixing is (u+v)/sqrt2, (u-v)/sqrt2.
The normalized analyses multiply both by another 1/sqrt2.
Thus the final per-copy entries are (u+v)/2 and (u-v)/2.
At reflected copy q' the u,v entries interchange; positive entries
agree and negative entries have opposite sign. At a fixed point u=v,
so the positive entry is v and the negative entry is zero.
This is exactly RETURN108-G's duplicate convention with multiplicity
represented by copies rather than a preweighted distinct-zero column.

For a nonfixed orbit of common multiplicity m, canonical weighted
coordinates are p=sqrt(m/2)(u+v), n=sqrt(m/2)(u-v).
The full normalized source has 2m copies per sign, equal to
p/sqrt(2m) for positive copies and +/- n/sqrt(2m) for negatives.
The embeddings with these constants are isometric, and disjoint orbit
supports give a global direct-sum isometry. At a fixed point there are
m equal positive copies v and zero negative copies; active positive
coordinate sqrt(m)v embeds by constant 1/sqrt(m).
No extra m or factor two is permissible.

## 3. Exact arbitrary finite packet weight

The actual selected projection accepts any finite set s of divisor
copies. Let k_i be the number of selected copies of a nonfixed orbit i
of common multiplicity m_i, and w_i=k_i/(2m_i) in [0,1].
Let n_i denote the canonical weighted negative sample above. Then

    ||selected normalized negative analysis||²
          = sum_i w_i |n_i|²,
    ||background normalized negative analysis||²
          = sum_i (1-w_i) |n_i|².

Fixed-point negative copies contribute zero. This follows literally
from the inspected coordinate projection and normalization; it does
not assume s is reflection-complete. The exact reduced finite
selected synthesis is M W^(1/2), where M is the canonical synthesis
for the finite set of orbits hit by s and W=diag(w_i)>0 there.
For each hit orbit the normalized signed constant vector on its
selected k_i copies is an isometric embedding of the one active
selected coordinate, with entries +/-1/sqrt(k_i). Thus the redundant
selected complement is again identity in the Picard marker, and the
active selected marker is

    C_s = (I+W^(1/2) M* Q^-1 M W^(1/2))^-1.

This is re-selection before inversion, not compression of the full
orbit marker. In the single-orbit case with full marker C=1/(1+k),

    C_s=C/[w+(1-w) C].

By contrast compression of the FULL duplicate marker onto the
normalized signed vector on the selected k_i copies yields
1-w+w C in that single-orbit model; the global limit is 1-w.
Compression after inversion can thus conceal a genuine crossing for
w<=1/2. Re-selection for any w>0 tends to zero if C does.
This generalizes 108-G from two duplicates to all analytic copies.
No actual historical compression error is claimed.

## 4. Actual source graph receives the witnesses without Green density

On the inspected full source graph, a compact smooth complex test
g on (-a,a) has finite log norm and hence a physical log coordinate.
The raw sample estimates of 108-C (uniform strip integration by parts
and multiplicity-weighted zero counting) give the two l2 hypotheses
of neutralActualZetaFullSourceLiftOfRawSamples, in the literal
ordinate convention z=-i(rho-1/2), since exp(i z x)=exp((rho-1/2)x).
Its physical coordinate is exactly g; no parity hypothesis or Green
density is needed. Its normalized sample entries are the isometric
duplicate embedding of canonical samples (choose negative orientation
consistent with u-v rather than v-u). Thus the conditional density
witnesses of 108-C enter this particular full source graph, at the
paper level, without the missing historical ConnesGreen predicates.

For each hit genuine orbit i in an arbitrary finite selected packet,
w_i>0. By global canonical density choose a smooth g with full
canonical samples arbitrarily close to e_i^-; the selected energy
approaches w_i and the positive energy approaches zero. A negative
selected direction therefore exists at a finite sufficiently large
aperture conditional on that genuine orbit. More quantitatively,
if full sample error <e, then positive norm <e and selected norm
>=sqrt(w_i)(1-e), by the contribution of i. Choose
0<e<sqrt(w_i)/(1+sqrt(w_i)) to get strict negativity. This bound
is independent of the graph norm. It does not locate a zero or give
an aperture. The full source quadratic identity attaches the full
signed form to the literal zero series and, for a smooth g, to its
C2 physical correlation, with the stated real-part convention.
No claim about arbitrary restricted parity carriers.

At the paper level, using 001G's bounded fixed-window raw analysis
on the corresponding log norm, every log vector has square-summable
samples and the full source graph is the graph of that bounded map;
its projection to the physical coordinate is a Banach-space
isomorphism at each aperture (constants may depend on aperture).
This does not identify the separate Hilbert graph norm with the
Dirichlet energy or prove Green synthesis surjectivity into that
graph. Avoid transplanting a Picard regularizer between these metrics
without checking inclusion and norm.

## 5. Disposition

The executable normalization and arbitrary selection operation have
now been recovered in the surviving IV091 full-source lineage, not
in the historical ConnesGreen actor implementation. We can use this
newer exact graph for independent witness transport without waiting
for the missing historical code. That is a model-specific attachment,
not evidence that the old code has identical columns or metric.
Historical source custody remains blocked; stop routine searches
unless a new archived ref is supplied. Do not conflate full graph
with the Green-generated closed graph subspace.

No new runner or formal theorem promotion. Next investigate actual
Green-generated range density in the physical Dirichlet/log norm
versus full graph topology, with explicit metric/source hypotheses,
or package the cross-route obstruction for independent review.
RH remains unresolved.
