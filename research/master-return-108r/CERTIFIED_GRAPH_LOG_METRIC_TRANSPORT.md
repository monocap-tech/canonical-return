# RETURN108-R — Certified graph receivers and exact cross-lineage log weight transport

Status: five graph-density declarations certified. Mellin bank
dependency build still running. Cross-lineage transport derived
analytically from literal source definitions; not a compiled bridge.
No global positivity or RH result.

## 1. Native receipt

Weil-Lab graph run 45926487392 job 140268875004 success.
Workflow commit 41c5c68e4cc0c1bbdf14902021772185781819fb.
Base GraphDensity blob 7ac84f8fa2ccb3dc746651d2cb420bb576ac8b1f;
corrected PacketGraphDensity blob 2b5554dd2dc5bd8a5114f286114033bec60f68c7
at ed7d208bc60f754784a4f8e14776f7d36e26119b.
Reopened logs: 3101 build jobs and five exact named PASS lines:
dense_image_of_dense_packets, closed_graph_eq_univ_of_dense_physical,
closed_graph_eq_univ_of_packet_transport,
graph_top_and_surjective_of_closed_packet_bound,
packet_bound_extends_to_closure. Final strict audit PASS 5/5.
Artifact 15167045951, return108q-lean, 1712 bytes, reported
SHA256 374e2ed51d6471c93dee11e64e25ead4e2966bc067ebbfd837e50a22f3f2f278.
Not independently downloaded or rehashed in this pass.
Cumulative inventory 19/27 native certified. Eight Mellin candidates
outstanding. Corrected Mellin run 45926491992 job 140268889502
has passed preflight, dependency resolution and cache; dependency
build still running at last read. No compiler conclusion on repaired
V10B5c or MellinKernelBank yet. No new runner activated.

## 2. Literal cross-lineage source audit

IV091 pin f2e8ec89ae4bf8687a652b299b730d6d1baa43f7:
Arithmetic/LogarithmicForm.lean ae932d4f3f7e4d6d06dbd778a5708927172238dc
defines w_I(xi)=log(exp(1)+|xi|).
Morphology/NeutralLogHilbertCarrier.lean
679d58b8f0cfa6b4f7d28c48ce844a41cdce2844:
neutralLogHilbertToCanonical, neutralCanonicalToLogHilbert and their
two inverse identities; exact neutralLogHilbertCarrier_norm_sq.
NeutralLogSourceDomainRealization.lean
e36c1e9c4cbc08ca4092701afaea3329873b738e:
neutralLogHilbertCanonicalLinearEquiv.
NeutralLogFormEnergy.lean 8b91dba883fc31e1e77645cc6e8bb4:
literal weighted Fourier mixed and norm-squared identities.
ActualZetaGreenCanonical.lean db60a9a6796e6b83cc066e92d41ff13d9f28e33e
supplies Green synthesis weak derivative and physical support.
ActualZetaGreenSampling.lean 3795b6fd172708c9b3a0549598429df1fcfdca58
defines actual columns, square summability and physical synthesis.

RECON repair head 8e0561ae553a924a9d91c68f38cd992ec607d7c3:
V08I logWeight is w_R(xi)=log(2+|2*pi*xi|).
V08J.physicalInclusion_injective, physicalInclusion_range_iff,
liftPhysical and physicalInclusion_lift give the exact supported
L2 carrier characterization and reconstruction.
Both use the Mathlib Fourier variable xi, not different independent
physical functions. The literal weights differ, including at xi=0.
The IV091 NeutralLogHilbertCarrier and ActualZetaGreenPacketDensity
files returned 404 on this RECON repair head. A native combined
import is not available merely by adding a file on this branch.

## 3. Explicit uniform weight and metric comparison

Write r=|xi|. Since e<3, 2*pi>1 and 2*pi<8, for r>=0
e+r <= 2(2+2*pi*r), and 2+2*pi*r <= 8(e+r).
Taking logs:
    w_I <= w_R + log 2 <= 2 w_R,
    w_R <= w_I + log 8 <= 4 w_I,
where w_R>=log 2 and w_I>=1, log 8<3.
Moreover |w_R-w_I| <= log 8. These bounds are independent of aperture.
For the same supported L2 vector f, finite weighted energy in either
norm implies finite energy in the other by measurable domination.
Both physical ranges therefore coincide. Define the unique
same-physical-vector map U_a from IV091 carrier X_I to RECON X_R
by V08J.liftPhysical on the IV091 physical vector using the bound.
Define inverse by neutralCanonicalToLogHilbert on the RECON physical
vector using the reverse bound. Literal physical reconstruction and
injectivity of the two maps imply they are mutual inverses and
complex linear. Exact norm-squared identities yield
    (1/2)||g||_I^2 <= ||U_a g||_R^2 <= 4||g||_I^2.
Thus ||U_a||<=2, ||U_a^-1||<=sqrt(2). No isometry is claimed.
In the ambient weighted Fourier coordinates the map is pointwise
multiplication by sqrt(w_R/w_I), with reciprocal inverse, as follows
from the actual reconstruction identities. Support is preserved in
physical space, not by an assumption on support in frequency space.

For the packet graph bound, if the canonical RECON unsigned orbit
sample map has norm B_a and its exact divisor-copy isometry is
attached, the IV091 unnormalized P/N graph energy is twice that
sample energy. The same physical vector therefore satisfies
    ||g||_G^2 <= (1+8 B_a^2)||P_I g||_I^2.
This is the correct conservative bound when B_a is expressed in
the RECON norm rather than IV091 norm. All constants finite for fixed
a. The graph receivers just certified use this bound on the packet
span and its dense physical image to conclude graph-top and onto.
The coefficient attachment is still a separate native obligation.

## 4. Compact metric defect

At the paper level the mixed norm difference on the shared physical
domain equals integral of (w_R-w_I) times conjugate(Ff)*Fg.
The multiplier delta=w_R-w_I is real bounded by log 8 in modulus.
Let L_I be physical inclusion from X_I into L2(R); on a fixed aperture
it is compact by the low-frequency Hilbert-Schmidt truncation and
logarithmic high-frequency tail estimate in 108-J. Let F denote
the unitary Plancherel map and M_delta its bounded multiplier.
Then
    U_a^* U_a - I = L_I^* F^* M_delta F L_I
is compact. This is an equality of operators on X_I from their
mixed forms, not an assertion that U_a itself is compact.
If the RECON signed covariance D_R=I+K_R with K_R compact and the
actual same-vector sample attachment identifies the signed forms,
then D_I=U_a^* D_R U_a=I+(U_a^*U_a-I)+U_a^*K_R U_a is also
identity plus compact. No uniform coercivity or nonnegativity
is inferred. In particular metric changes do not obstruct this
full-kernel completeness route. Actual signed covariance
attachment and native compact-defect proof remain separate.

## 5. Disposition

The generic graph-top transport is certified, the exact metric
bridge is now specified at the analytic level. Do not conflate
the branches or assume their definitions definitionally equal.
Next: read the corrected Mellin build. If the inherited dependency
still blocks elaboration, split the finite observation kernel
receiver away from the budget imports, retaining the literal
coordinate hypotheses. Then plan a combined pinned IV091/RECON
source integration or prove the same interface directly in IV091.
No new native declaration or workflow activated in this pass.
RH unresolved.
