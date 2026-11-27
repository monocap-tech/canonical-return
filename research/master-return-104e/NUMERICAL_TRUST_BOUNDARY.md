# RETURN104-E — Numerical trust boundary and proof-package acceptance ledger

Status: source-level audit and acceptance plan; no new execution or theorem
promotion. RETURN104-C run 45917729699 is still executing. RETURN104-D is
staged at 267a3b32d81070fa0303c146753f0ef92e87d794, not executed.

## 1. Which arithmetic is actually still a premise?

The reconstructed high-band (RETURN089/090) and odd strict-tail auxiliary
Gram (RETURN088) use Python Fraction and NumPy object integer endpoints.
Their interval matrix product computes all four endpoint products for each
inner index, sums min/max separately, then divides with directed integer
rounding. No binary64 dot product is in the generating arithmetic of
these three scripts. The historical return073/upgrade routines imported
for reference containment are not the generating arithmetic.

The full and restricted finite high bands consume the original
RETURN086 direct_sums.json at all powers 2 through 46 and the four
linear features 1,log x,S,C. The ten symmetric scalar product indices
are [1,L,L²,S,C,LS,LC,S²,SC,C²], in the literal order encoded in the
col matrix of both scripts. The 101-A Arb direct sum regenerates
exactly those families with the identical K=32768, U=262144 and
strict split k>65536, and 103-B binds its saved enclosures to the
literal historical scalar bank. It is not necessary to certify
historical Abel recurrence arithmetic in addition to the direct
sums in order to use these reconstructed matrices.

The remaining D0 import is ref.D0. Reopened the actual archived
return073_band.py (not inferred from the import name). Its definition is
exactly the Q=10^55 directed src.endpoints(gR+src.G0+4), with
r=2*iv.pi*K/5 and gR=src.G0+fast.sums((0,K))[0]+4*r²/(r²+1/4).
This is the same endpoint expression that RETURN101-A compares
against Arb at K=32768. The two high-band scripts use the explicit
rational pi interval [314159265358979323846,314159265358979323847]/10^20.
This is wider than the known exact pi; document and independently
assert containment rather than assume equality to RETURN062's interval.
The actual sin/cos polynomial bounds |S|,|C|<9 and log(8)<3 justify
the positive remainder used in cross/restore. For strict auxiliary
tail the trapezoid correction uses -1/2, not +1/2.

These findings distinguish reference-only imports from mathematical
inputs. They do not claim that every imported historical routine
can be erased from the reproduction environment. In particular the
source-level D0 identity above should receive a literal automated
assertion in the final proof-package replay.

## 2. Exact integer backends

RETURN104-B successfully replaced the finite auxiliary Gram midpoint
products of RETURN092-B/C by fmpz_mat multiplication. All 815414
entries of the common even, common odd and augmented odd matrices
are exactly equal, and all four full terminal inputs are unchanged.
The exact backend uses the same interval radii and tail remainders.
The source-dependent matrix identity is not established by checking
integer-valued binary64 output alone.

RETURN104-C similarly targets the RETURN087 finite and compressed
feature dot routine. Until its exact bank equality succeeds, retain
the IEEE binary64 exact product/accumulation assumption for these
banks; do not mark it removed. Confirm the actual archived dot
function and the feature driver's imported helper, not just the
name selected by the API probe.

RETURN104-D is the remaining direct substitution of the actual
infinite leading tail phase chords, weights and odd Gamma anchor.
RETURN101-A's separate 128-key distance test is not by itself proof
that the particular chord lower endpoints used by RETURN090-B and
091-C are valid. Exact equality of the re-evaluated rational banks
would establish the needed transfer; if equality fails, examine the
specific intervals and either prove enclosure of the original bank
or replay the downstream terminal certificate using independently
valid replacements. Do not treat overlapping intervals as a transfer.

## 3. Minimal reproduction trust model

For a conventional computer-assisted fixed-aperture proof the
remaining software TCB should be stated explicitly: Python exact
integers and Fraction, exact object-array operations, FLINT integer
matrix arithmetic, and Arb rigorous real/complex balls. Python-flint
0.8.0 and NumPy 2.2.6 are pinned in recent replays. Original
reference-reader dependencies (mpmath, SciPy, NumPy) should not
be counted as mathematical premises when the actual consumed
interval is independently validated and the reference check is
only a comparison. Hashes authenticate byte identity, not numerical
soundness; CI success certifies execution of the assertions, not
a proof assistant theorem. Document the exact build, compiler and
runtime assumptions rather than claim zero software trust.

The analytic TCB includes the standard Fourier/Parseval and digamma
identities and the Guinand–Weil explicit formula for the spectral
corollary. The physical inequality does not require RH or that
corollary. RETURN103-A/C provides an assembled implication and
internal review, not independent referee acceptance.

## 4. Acceptance sequence

1. Confirm RETURN104-C output equality and exact backend coverage.
2. Execute staged RETURN104-D only after the current runner is done;
   inspect any interval differences without relaxing the terminal targets.
3. Preserve successful RETURN103-B and RETURN104-A/B/C/D archives
   with digests and immutable-by-policy repository custody. Current
   known accepted supplemental artifacts: 15162640511 (103-B),
   15162478276 (104-A), 15161989586 (104-B).
4. Bind new exact backends to original consumed banks; verify literal
   high-band D0 and pi endpoints; document the complete arithmetic TCB.
5. Consolidate analytic derivations and dependency ledger into a
   reproducible fixed-aperture proof package. Keep the source
   statement and limitations explicit. Independent mathematical
   review remains a separate acceptance standard; native Lean
   formalization remains a separate workstream.

Stop further redundant audits once the stated dependencies have
been resolved, except where a specific new flaw is identified.
No increase in aperture or support globalization follows from
removing a software dependency. Conditional floor 3*10^-59 unchanged.
