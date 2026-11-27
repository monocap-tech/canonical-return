# RETURN106-A — Mathematical sufficiency review and final attachment criteria

This is an internal mathematical review of the integrated RETURN105-B
package (42a7ec0b649e2301cb9f25dc5581a292b29e8a7e), not an external
referee report. No numerical optimization, theorem promotion or new
runner. The permanent supplemental custody commit is
4399c20b2b708260463dffa25f02797183d1983f; manifest Git blob
47e3d9d70fdda55f4b0b2c76b1384101e394915d and high-band attachment
SHA-256 7b3323ec05f9f476c33150d3b3c400ef69f0ffdc8c63c5cca78ab5b41b78a21c.

## 1. Receiving implication: sufficient hypotheses

For one parity let H be the actual retained form in its rational
coordinates, D_phys=(5/4)D the complementary comparator restricted
to ker C0, and t the actual mixed functional. The following hypotheses
are sufficient and do not refer to RH or any zero-location assumption:

(i) H >= h I with h>0, and physical retained mass <= M ||x||²;
(ii) D_phys >= d ||q||² on the actual complementary kernel,
     d=297/6250;
(iii) the modeled source m has squared dual norm <= mu H[x];
(iv) the exact source difference t-m has dual norm <= sum eps_i sqrt(H[x]);
(v) f=p+q is an orthogonal physical decomposition in the form domain,
    and W(f)>=H[x]+D_phys[q]+2<t(x),q>.

Proof: define the dual norm variationally as
sup_{q != 0} |<v,q>|/sqrt(D_phys[q]). By coercivity it is finite for
the bounded l2 functional v and obeys triangle. Thus t has dual norm
<= r sqrt(H[x]), r=sqrtcap+sum eps_i. The elementary inequality

    H+D_phys-2r sqrt(H D_phys)
     = (1-r)(H+D_phys) + r(sqrt(H)-sqrt(D_phys))²

gives (1-r)min(h/M,d)||f||², including the cases p=0 or q=0.
The coefficient r is nonnegative and <1. Orthogonal parity splitting
gives the common minimum. This is the complete receiving proof; there
is no requirement that the saddle K itself be positive, or that the
three error bands be orthogonal in the inverse metric.

## 2. Independent analytic checks in this review

Reopened RETURN098-A/B/C, 099-A, 102-A, 090-B, 091-C, 091-B,
092-A and the literal terminal checker. The physical exponential
bilinear formula uses complex bilinear, not sesquilinear, products.
The m=0 triangle takes its continuous divided-difference limit;
the m>=1 denominators are nonresonant for real exponent parts
between -1/2 and 1/2. The retained source evaluator 101-C uses
the corresponding digamma divided differences and explicitly pays
the trigamma and exponential decay remainders.

The drift estimate is valid: the derivative kernel
4r beta/(beta²+r²)^2 has integral 2/r and total variation
9/(2 sqrt(3)r²). For sample beta=2m+1/2 in cells of length 2,
the sum minus half-density integral is bounded by this variation.
Consequently |g'(r)-1/r|<3/r², so the weaker 8/R used in the
actual leading-tail contraction is safe.

For h=sum_{m>=1} (beta_m²+r²)^-1 the missing lower integral
relative to pi/(4r) is <=5/(4r²). The positive geometric
bounds on r² e and r² f are <1. Thus the quoted 6/r alpha
error and 8/R+3/r² d error follow. For the finite bands
|alpha0|<11, |d0|<23 because |D0|<10, 2|Pc|<9 and log x<3.
All coefficients of the six error envelopes are nonnegative.
The 102-A error script uses full charted norms and a geometric
remainder bounded by the sum of raw-row norms times coefficient
absolute values divided by 1-rho. It expands the square of
the infinite envelope termwise, including both off-diagonal
copies. Its power-log integral reserve 2q! is justified by
the total variation of x^-p log^q x on [1,infty), q=0,1,2.
No use of the older under-rounded scalar norm cap is needed.

The odd infinite model is
l B0 - h_s G', where l=(2Pc-g(r)+g(R))/r and
G'=5(D0 B0-2L0)/(2pi). RETURN090-B actually evaluates the
cross Gram with the negative sign; it pays drift in l² and
l-crosses, then includes the strengthened restoration cross.
The compensation error uses
|r/(r²+1/4)-1/r| <= 1/(4r³).
No omitted log drift in the exact source error is found.
The exact finite and leading model supports meet at the cutoffs
without double counting. The restricted band is a weight change,
not extra coordinates. RETURN091-B applies lambda=174/1199 to
the tail TT/U and separately reads the strengthened restoration
cross; RETURN092-A consumes its exact serialized TT/U. Its
workflow pins the 091-B ZIP digest, and 103-B preserves the
resulting source fields through the auxiliary/retained replacements.

The restricted inverse comparison is variational. Since
D>=A-Gneg*Gneg on the kernel, extending the supremum to
all l2 and dropping the positive updates yields the upper dual
bound. Woodbury and the trace bound on Gneg P A^-1/2 yield
the stated restricted factor, without assuming commutation.
The inverse Neumann cap and trial residual in the terminal
checker bound the indefinite saddle error by its absolute
operator norm. The weighted retained radius majorant acts
before the ill-conditioned chart congruence. These arguments
support (i)-(v) once actual interval membership is attached.

## 3. Specific final evidence interface, not a new search

The first 400 mixed rows are a distinct interface. RETURN104-A
strictly encloses them inside the authenticated RETURN097-B
alternate bank. RETURN097-B compared against return067_contract
rawrows(par,S); the independent RETURN087 finite contraction
used its own reader for the initial 24 normalized and remaining
376 raw/chart rows (see original MASTER_RETURN_087.md).
104-C proves equality of the 087 OUTPUT contraction banks, not
identity of those two INPUT readers. Later grouped rows in
104-A compare directly to the original archived files, and
the compressed features compare directly to return080/features.json.

The accepted source lineage supports the intended connection,
but for a standalone literal proof package require a short
deterministic check that the 104-A early interval bank lies
inside the actual 087 finite input rows at the identical
scale, labels, chart and rephasing. Alternatively identify
an already executed assertion establishing that reader identity
in the authenticated 087 source. This is a targeted attachment
criterion; mere overlap of the two readers is insufficient.
It is not a discovered counterexample or a failure of the
physical inequality. If a difference is found, use the actual
Arb early inputs and replay their finite contribution and final
terminal comparison without changing the threshold.

The 103-B joint binding has already checked the terminal
field-preservation chain and physical floors. The 104-B/C/D
exact replacements and 105-A attachment remove the identified
numerical backend issues. No further routine independent
re-generation should be mandated once this input interface is
documented and the evidence index pinned. The proof still
requires ordinary mathematical scrutiny of the analytic lemmas
and software soundness, as with any computer-assisted proof.

## 4. Globalization is a different mathematical obligation

The fixed-aperture statement bounds all real tests in any
translated interval of length 5/2, subject to acceptance.
It does not follow that sums of such tests on many intervals
are positive. For a decomposition f=sum f_i, local inequalities
control the diagonal W(f_i) but leave the cross W(f_i,f_j).
Even a two-dimensional symmetric Gram with diagonal c>0 and
off-diagonal 2c is positive on each coordinate axis and negative
on their difference. This is a generic logical obstruction,
not a counterexample to the actual Weil form.

The next genuinely global object must control the off-diagonal
translated-support interaction (or provide a different exact
global transport). The prime-shift set and pole contribution
change as the support diameter grows. In particular a partition
of unity or rescaling argument without cross-term control is not
a globalization proof. Do not infer RH from the local certificates.

Acceptance state: the abstract receiving theorem is valid;
the reviewed analytic derivations show no detected contradiction.
This is internal review only. The explicitly identified early
input interface, a final evidence index, and independent
mathematical review remain. Conditional floor 3/10^59 unchanged.
