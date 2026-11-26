# RETURN098-B — Tail compensation and restricted inverse receiving lemma

Scope: analytic proof of the transport identities for the fixed RETURN083
model. This is not a fresh interval certificate, native Lean proof, or
unconditional sign theorem.

## 1. Three disjoint source-error supports

Let J be the complementary even grid of RETURN096-A. Set
K1=16384, K2=32768, K3=262144. The original mixed source t is the
sequence of rephased pairings of RETURN098-A, in the retained chart.

The model m has four pieces: exact rows k<=K1; 36-feature model on
K1<k<=K2; 24-feature model on K2<k<=K3; and leading asymptotic model
on k>K3, with the odd compensation described below. In particular
k=32768 belongs to the variable band and uses w=150000 (A=87/500);
k=65536 uses the old A=87/500; k=262144 uses A=519/1000 in odd
parity. The cutoff inequalities of RETURN096-B justify these choices
without assigning any coordinate two weights.

Let P_i denote coordinate projections onto the three nonexact pieces.
Then t-m=sum_i P_i(t-m). There is no assertion of orthogonality in the
constrained inverse norm. The triangle inequality in that very norm,
not the sum of squared errors, is the receiving step.

## 2. Odd gauge is compensation, not a null functional

Write r_k=2*pi*k/5, c=25/(16*pi^2), h_s(k)=k/(k^2+c).
Let B0 and L0 be the literal leading retained vectors from the
normalized raw coefficients (L0 includes c_raw+d_raw/2 in the
notation of RETURN066). Define D0=g(R)+g(0)+4 at R=r_K3 and
G=D0 B0-2 L0. Then G'=5 G/(2*pi) and

    r_k G/(r_k^2+1/4)=h_s(k) G'.

On the entire complementary grid h_s is a constraint row, but
P_3 h_s need not be a constraint row. Hence one must not discard
the tail-only compensation by appealing to sinh orthogonality.

With H_3=P_3 h_s, the odd leading source on k>K3 is the
rank-one sequence l(k) B0, l(k)=[2 P_c(r_k)-g(r_k)+g(R)]/r_k.
The literal odd tail model is

    m_3 = l B0^T - H_3 (G')^T.

This is the model in RETURN070's leading construction, including
the negative cross terms. The actual tail error is

    t_3 - m_3 = (t_3 + H_3 (G')^T) - l B0^T.

Thus the compensated error estimate is an estimate for the
difference to the actual stored model. No constraint term has
been removed from just one side of an inequality.

Let m_f be the combined model through K3 and let W=A^-1 be the
diagonal reciprocal of the final piecewise comparator. Its supports
are disjoint from l and H_3. With scalar d=<l,Wl>,
e=<H_3,WH_3>, a=<l,WH_3>, the complete Gram is exactly

    TT = m_f^T W m_f + d B0 B0^T + e G' G'^T
                         - a (B0 G'^T + G' B0^T).

The complete auxiliary cross is

    U = T W m_f + (T W l) B0^T - (T W H_3) G'^T.

These identities include the restoration rows of T, whether or not
a row is a constraint. They hold unchanged for a strengthened
piecewise A, provided all scalar products use its actual reciprocal.
They explain why the cross between the two asymptotic pieces must
be kept despite their being stored in different historical banks.
RETURN070 source and the RETURN068 gauge implementation are
consistent with this sign convention; the literal K3 gauge receipt
and the later signed reweighting still require final joint audit.

## 3. Restricted signed inverse error bound

Work in the real l2 grid with diagonal A bounded above and below by
strictly positive constants, as in RETURN083. All constraint and
signed Cauchy rows are bounded functionals. Let
D_-=A-Gneg^*Gneg, where Gneg uses the normalized negative rows.
RETURN095-B certifies ||Gneg A^-1 Gneg^*|| <= rho=13/25 for the
baseline A. Increasing A can only decrease this Gram. Thus D_-
is coercive and

    D_-^-1 = W + W Gneg^* (I-Gneg W Gneg^*)^-1 Gneg W.

The actual receiving operator S is the inverse of the positive
complementary comparator on the constraint kernel, extended by
zero on the orthogonal complement in the variational pairing.
Positive signed updates, positive restoration, and restriction of
the variational supremum give, for every source x,

    <x,S x> <= <x,D_-^-1 x>.

This follows directly from sup_q [2<x,q>-D[q]], and requires no
commutation between coordinate projections and the inverse.
For x supported in a band E, put P=P_E and
tau_E >= tr(Gneg P W P Gneg^*). Then

    ||Gneg W x||^2
       <= ||Gneg P W^(1/2)||^2 ||W^(1/2)x||^2
       <= tau_E <x,W x>.

The inverse of I-Gneg W Gneg^* has norm <=1/(1-rho).
Consequently the precise restricted error receiving inequality is

    <x,S x> <= (1+tau_E/(1-rho)) <x,W x>.

The unrestricted version is <x,S x><=(1-rho)^-1<x,W x>.
Both are valid for x of either sign, and for the difference of
two sources. The physical inverse multiplies S by 4/5.
Dividing by the positive H chart floor gives exactly the form
of the bandwise squared error estimates in RETURN085. It remains
necessary to check the numeric tau_E and source norm inputs in
the identical chart and with the final A; this derivation alone
does not establish those particular numbers.

## 4. Assembly boundary

These are exact analytic identities with the following dependencies:
the actual retained chart, the exact source-to-model piecewise
definition, the normalization of the signed and restoration rows,
the reconstructed TT/U/K intervals, and the certified error reserves.
Use RETURN095-A only once those dependencies are attached.
Do not promote an overlapping interval as contained, treat the
tail gauge as zero on the constraint kernel, or replace a bandwise
inverse estimate by an unsupported orthogonality claim.

No runner is activated by this note. RETURN097-A/B/C remain the
finite scalar regeneration gates. No native Lean, arbitrary-support
positivity, or RH conclusion.
