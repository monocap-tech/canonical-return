# RETURN108-I — Green completeness in Dirichlet energy and full source graph

Status: internal analytic density theorem using classical unconditional
zero counting and distinct critical-line zero supply. Not a native Lean
verification or historical ConnesGreen source recovery. No RH claim.

## 1. Source attachment

Reopened Weil-Lab IV091 pin f2e8ec89ae4bf8687a652b299b730d6d1baa43f7:
WeilDefect/Arithmetic/ActualZetaGreenPacketDensity.lean
2325cfc5812ee50a993de4e4a382b7cc0e3d532f,
ActualZetaHilbertSourceGraph.lean
595316a5c06f79746a8ac838cef1463d029fc3df,
ActualZetaGreenClosedSubspace.lean
fca4b0cbe82c109976b2692dbe5d6e60c6309f1c,
ActualZetaGreenGraphClosure.lean
0390929f1294b522f82b552ecb7fa0d9ba76f909,
ActualZetaGreenSampling.lean
3795b6fd172708c9b3a0549598429df1fcfdca58,
WeilDefect/Morphology/NeutralLogSourceDomainRealization.lean
e36c1e9c4cbc08ca4092701afaea3329873b738e.
Use literal index, graph and normalization sources pinned in 108-H.
Existing PacketClosedSpan_eq identifies packet closure with Green graph
closure, NOT the entire graph. The source graph norm is the sum of the
physical logarithmic norm squared and both UNNORMALIZED P/N norms squared.
The existing graph-closure mixed and diagonal form attachments apply on
that closure, not automatically to all source vectors. These proof bodies
were inspected; no fresh compilation/axiom audit.

## 2. Dirichlet completeness by exponential-type uniqueness

Fix a>0. Let V=H0^1(-a,a;C) with energy inner product
int conjugate(f')g' + conjugate(f)g/4.
The actual IV091 column indexed by rho solves
L G_rho=exp(-i z_rho x)=exp(-(rho-1/2)x), with L=-d²+1/4,
z_rho=-i(rho-1/2), and zero endpoints. This is the literal
DirichletProblemOneColumn expression, verified mathematically by
applying L to its particular solution and homogeneous boundary terms.
The denominator is nonzero on the open strip. Equivalently relabel
the full family by mirror rho'=1-conj(rho) to obtain source
exp((conj(rho')-1/2)x), as in 108-E. The relabeling is bijective
by zero reflection symmetry and is not a per-index column equality.
Each column lies in V. Multiplicity copies do not affect the span.

If h in V is energy-orthogonal to all columns, weak integration
by parts (first on smooth compact tests, then by H1 density) gives
0=<G_rho,h>_D=int exp(-(conj(rho)-1/2)x) h(x) dx.
By the bijection rho->1-conj(rho), its entire transform
H(z)=int exp(zx)h(x)dx vanishes at all rho-1/2.
The integral is entire by dominated differentiation, and
|H(z)|<=||h||_1 exp(a |z|).

If H is not identically zero, choose z0 with H(z0)!=0.
Jensen on the disk about z0 of radius 2R gives
n_(z0)(R) log 2 <= log max_|z-z0|=2R |H(z)| - log |H(z0)|
<= 2aR + O_(h,a,z0)(1). Therefore the number of distinct zeros
in |z|<=R is O(R) as well (use radius R+|z0| around z0).

External unconditional input: Conrey, "More than two fifths of
the zeros of the Riemann zeta function are on the critical line",
J. reine angew. Math. 399 (1989), 1–26, Theorem 1 p.4,
https://aimath.org/~kaur/publications/24.pdf .
The original printed theorem was checked. It gives a positive
liminf proportion of simple critical-line zeros (in fact >=0.401).
Together with Riemann–von Mangoldt, this gives distinct simple
critical zeros up to height R numbering >> R log R. This is the
needed distinct-zero supply, not a count with duplicates treated
as extra zeros of H. Even just these critical zeros contradict
Jensen's O(R). Thus H identically zero and h=0 by Fourier
uniqueness on L2. The energy orthogonal complement is zero:
the finite complex span of the actual Green columns is dense in V.
In fact simple critical-line columns alone suffice. No RH or
lower separation bound, no formalization of Conrey is claimed.
This is the opposite direction of density to 108-C (there:
physical tests dense in active coefficient space as support grows).

## 3. Fixed-window bounded raw analysis, direct analytic proof

Let X be the supported logarithmic carrier, using angular Fourier
variable with weight w(y)=log(2+|y|). In Mathlib frequency xi
this corresponds to y=2 pi xi; if the original logarithmic
weight uses log(2+|xi|), it is equivalent by constant factors.
Use the actual original carrier via the inspected canonical
physical/log identification. No equality of norms for differently
scaled Fourier coordinates is presumed.

For f in X its centered bilateral transform F(z) is entire
since f supported in [-a,a] and L2, hence L1. Each centered zero
a_rho lies in |Re z|<=1/2. The analytic submean inequality on a
radius 1/2 disk gives
|F(a_rho)|² <= (4/pi) int_disk |F(z)|² dA(z).
Summing with analytic multiplicities by Tonelli, all disks lie
in |Re z|<=1; at a fixed ordinate y their centers lie within
distance 1/2 of y. The unconditional unit-height zero count
N(T+1)-N(T)=O(log(2+T)) (and conjugation for negative heights)
bounds the number covering any point by C log(2+|y|), including
bounded low heights. Consequently
sum_rho m_rho |F(a_rho)|² <= C' int_-1^1 int_R
                                  w(y)|F(b+iy)|² dy db.

For each b in [-1,1], choose a fixed smooth compact cutoff chi
equal to one on [-a,a], and put chi_b(x)=chi(x) exp(bx).
Then F(b+iy) is (up to Fourier sign) the Fourier transform
of chi_b f. The weight satisfies
w(y)<=w(y-v)+w(v)<= (2/log 2) w(y-v)w(v).
The first inequality follows from 2+|y| <=
(2+|y-v|)(2+|v|), the second from both weights >=log 2.
Thus weighted Young's inequality in L2 yields
||sqrt(w) Fourier(chi_b f)||_2 <= C
 ||sqrt(w) Fourier(chi_b)||_1 ||sqrt(w) Fourier(f)||_2,
with the Fourier normalization constant included in C.
The first norm on the right is uniformly finite in b:
chi_b has fixed compact support and all fixed-order derivatives
uniformly bounded, hence uniform rapid Fourier decay.
The same proof by L2 approximation justifies the convolution
for non-smooth f. Integrating over b proves the bound
sum_rho m_rho |F(a_rho)|² <= B_a² ||f||_log²
for a finite aperture-dependent B_a.
This uses the classical unit-height counting theorem already
pinned in historical IMPORTED_SOURCE_PINS EXT-3 (Titchmarsh,
2nd ed., theorem 9.2). It is a standalone analytic receiving
estimate, not a newly compiled theorem or uniform-in-a bound.

The reflected sample family obeys the same bound by reindexing
the divisor involution. By the literal first 1/sqrt2 mixing
the squared unnormalized P/N norms sum to twice the raw norm²:
each full divisor copy is retained. Therefore the Hilbert source
graph lift L_a has norm bounds
||f||_log <= ||L_a f||_graph <= sqrt(1+2 B_a²)||f||_log.
It is a bijection because the graph equations uniquely determine
both source coordinates, and 108-H's actual lift accepts both
raw square-summability hypotheses. The constants absorb harmless
equivalence of weight convention. No norm equality or uniform
global bound is claimed.

## 4. Transfer of Green density to the exact graph

108-E proves that the zero-extended Dirichlet carrier maps
continuously and densely into X: Fourier log weight bounded
by constant times 1/4+y², and inward dilation/mollification
proves smooth compact tests dense in X. Section 2 gives
Dirichlet norm approximation of each smooth compact test by
finite sums of the original Green columns. By Section 3
their log convergence implies Hilbert graph convergence.
Each such finite sum is the actual finite-packet Green lift:
its physical coordinate is the finite sum of literal Green
columns, and both graph sample coordinates are uniquely
determined by that physical vector. No redefinition of a
column or of the graph norm is used.
Therefore every smooth compact graph lift is in the Green
graph closure. Since these lifts are dense in the entire
source graph by the log approximation and bounded lift,
we obtain at the paper level

    closure(span actual Green packets) = full Hilbert source graph

for each a>0. In IV091 terminology the target is
neutralActualZetaGreenClosedSubmodule a ha = top.
This is a mathematical statement about the interpreted literal
source objects, not an asserted Lean declaration of this equality.
In particular no surjectivity of the compact l2 Green synthesis
itself follows, nor any uniform control of synthesizing coefficients.

Combined with existing GreenGraphClosure_diagonal/mixed, this
would extend the original attached native Weil quadratic identity
to every vector of the full source graph on each fixed aperture.
This is a transport completion; it makes no assertion of
nonnegativity on the Green packets or full graph. If a genuine
off-critical negative channel exists, 108-H's smooth finite-aperture
negative witness is now approximable by finite Green packets
in the full graph norm on its containing aperture. Since the
bounded signed quadratic is continuous, sufficiently accurate
finite Green packets also have negative signed form. No actual
negative witness or coefficient choice has been identified.

## 5. Boundaries and next task

The only external arithmetic inputs here are classical distinct
critical-line supply and zero counting, not RH-strength bounds.
The new completeness argument should be independently reviewed
and its formal dependencies made explicit, especially the
Conrey/zero-counting import (or an already formalized weaker
superlinear distinct supply), Jensen uniqueness, and the
weighted sampling bound. Do not count the source-level target
as native Lean PASS. The original historical ConnesGreen code
remains unavailable, and this proof does not identify its metric
with the IV091 graph norm. The fixed-aperture RETURN numerical
candidate and external review are unchanged. No runner activated.

Next: audit whether the current certified RECON sampling and
critical-supply interfaces can discharge these dependencies,
and package the exact native graph-top target without introducing
an axiom or silently importing an unproved distinct-zero result.
RH remains unresolved.
