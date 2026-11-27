# RETURN103-A — Fixed-aperture receiving theorem and dependency ledger

Status: assembled computer-assisted theorem candidate, with explicit
certificate premises. This is not a native Lean theorem or an RH claim.
No arbitrary-support conclusion is available. The permanent evidence
commits are 87c9cd17b29da9d0a57cd426e0aa76a24eea0f65 (RETURN088–099)
and bce3d5491d0cabe684e05237ca3e0533f22a5914 (RETURN101–102).

## 1. Exact statement and domain

Let I=(-5/4,5/4). For real f supported in its closure put
F(z)=integral f(x) exp(i z x) dx, and c(t)=integral f(x)f(x+t) dx.
On the logarithmic form domain of zero-extended functions define

    W(f) = (1/(2*pi)) integral_R |F(r)|²
                  [Re psi(1/4+i*r/2)-log(pi)] dr
           -2 sum_n Lambda(n)/sqrt(n) c(log n)
           +2 F(i/2) F(-i/2).

The candidate conclusion is W(f)>=3*10^-59 ||f||_2² on this
domain. In particular for real C_c^infty(I), the unconditional
Guinand–Weil explicit formula identifies W with the sum of
F(z_rho) F(-z_rho) over the nontrivial zeta zeros, counted with
multiplicity, z_rho=(rho-1/2)/i. It is not |F(z_rho)|² away
from the critical line. This receiving theorem does not assume RH.
The explicit formula is needed for that spectral corollary, not
for the physical inequality itself. Translation of a real test
leaves W invariant (the pole moment factors cancel), hence any
translated interval of the same length is covered; dilation does
not preserve the prime translations or the Gamma multiplier.

## 2. Certificate premises

For each parity use precisely the authenticated 87-generator retained
space P and the rational chart C/T0 of RETURN059, as identified in
RETURN096-A. Let J be the complementary positive even grid with
16302 finite coordinates through 32766 and all even k>=32768.
Write the orthogonal decomposition f=p+q, p in P, q in P-perp,
with p represented by x in R^87 and q by rephased coefficients a
in the actual constraint kernel ker C0. Here C0 denotes the
constraint operator, not the retained chart C.

The following finite/infinite interval statements are the certificate
premises. Their actual enclosures and matrices must have the same
source/chart/coordinate labels, not merely overlapping numbers.

(A) H is the literal physical retained energy in the chart and
H[x]>=h||x||², with h the rational RETURN059 margin. The raw
matrix enclosure is the one used in RETURN093-B. The chart is
nonsingular, and ||p||²<=M||x||² with
M=(5/2)*87*||C/T0||_F².

(B) A is the literal piecewise diagonal of RETURN096-B/100-C.
On the constraint kernel, W(q)>=(5/4)D[a], where
D=A+L*J_s^-1 L includes the signed carriers and, in odd parity,
the eight restoration rows with metric R. The normalized negative
row operator G_- has ||G_- A^-1 G_-*||<=13/25. In particular
D>= (12/25)A>= (12/25)*(99/1000)I. A is bounded above as well
as below: the finite ledger is finite and the infinite tail is
piecewise constant. The physical energy D_phys=(5/4)D satisfies
D_phys[a]>=(12/25)*(99/1000)||q||² by Parseval.

(C) The literal mixed source t(x) satisfies W(p,q)=<t(x),a>.
Its model m(x) has the four disjoint frequency pieces specified
in RETURN098-B/099-A. For T_aux=(C0,L),
K=diag(0,J_s)+T_aux A^-1 T_aux*, TT=m* A^-1 m,
U=T_aux A^-1 m, the terminal certified matrices enclose these
operators in exactly these charts. K is invertible, and

    (4/5)(TT-U* K^-1 U) <= mu H,

with mu=49/50 even, 19/20 odd. The three errors e_i are the
literal restrictions of t-m, with dual norm relative to H[x]
bounded by the three RETURN102-A norm caps in physical units.
The tail-only odd compensation must be included in both model
and error, not discarded by invoking the global sinh constraint.

(D) All finite row, feature, moment, scalar and auxiliary banks used
in (A)–(C) are valid enclosures of their indicated mathematical
quantities. Their arithmetic implementation is directed interval
arithmetic or exact integer/rational calculation with the specified
rounding and remainder reserves. Neither downstream containment
alone nor overlapping independent intervals proves this premise.
Use the independently regenerated banks wherever they have replaced
the original, and require literal lineage equality at the terminal.

## 3. Receiving proof

The real physical form is symmetric and reflection invariant, so
split f by parity. Work in one parity. The retained modes and the
orthogonal complement belong to the log form domain whenever f
does, by the piecewise-smooth Fourier decay of retained modes.
The form is well-defined there by the bounded image, pole and
prime corrections to the circle log form. The mixed map extends
in the q variable to L2 by the square-summable log-weighted
Fourier coefficients of each retained mode. These claims are
derived in RETURN096-A/B and RETURN098-A.

By polarization and (A),(B),(C),

    W(f) >= H[x]+D_phys[a]+2<t(x),a>.

The constrained inverse formula of RETURN095-A applies: D is
coercive, K invertible, all rows bounded on l2. The saddle
response is A^-1(m-T_aux* K^-1 Ux); it lies in ker C0 and solves
the constrained variational equation. Thus its physical dual
square is (4/5)x*(TT-U*K^-1 U)x, bounded by mu H[x].
For each error band the variational comparison with A-G_-*G_-
and Woodbury give the unrestricted factor 25/12 or restricted
factor 1+tau/(12/25), multiplied by 4/5. This is the exact
receiving norm used in RETURN102-A.

By the triangle inequality, not an orthogonality assertion,
the full source dual norm is at most r sqrt(H[x]), where
r=sqrtcap+sum_i e_i, sqrtcap=99/100 even, 39/40 odd.
These rational caps square above the respective mu values.
The two certified r are <1. By Cauchy and 2uv<=u²+v²,

    W(f) >= (1-r)(H[x]+D_phys[a])
         >= (1-r) min(h/M,(12/25)*(99/1000))
              (||p||²+||q||²).

The sum of masses equals ||f||² by the orthogonal decomposition.
The same applies to the two orthogonal parities. RETURN102-A's
exact rational evaluation gives floors exceeding 3*10^-59 for
both. This proves the receiving implication on the specified
domain. Density of real C_c^infty in the supported log form
domain follows by contraction of the support and mollification
as in RETURN098-A; it does not enlarge the support length.

## 4. Dependency ledger

- Original physical identity and domain: RETURN096-A/B and
  RETURN098-A/C, analytic derivations; external explicit formula
  only for the spectral interpretation.
- Exact complement partition and constraints: RETURN096-A native
  Python integer/rational audit (not native Lean); 5/6 constraint rows.
- Signed symbol and restoration: RETURN096-B, RETURN094-A,
  RETURN095-B; all checks passed, analytic form normalization 5/4.
- Retained physical matrix and chart floor: RETURN093-B, RETURN099-B,
  RETURN101-C; independent Arb containment against eight banks.
- Exact finite mixed rows: RETURN097-A/B, 8111 per parity, and
  RETURN087 independent finite contraction; all containment checks passed.
- Variable-band model: RETURN097-C, RETURN087 feature contraction,
  RETURN102-B actual moment intervals; K=16384 obtained by exact
  4^n rescaling of the 32768-anchor moments. 8192 coordinates.
- High-band model: RETURN089/090 and RETURN101-A directed scalar
  sums and phase constants; 114688 full, 98304 restricted coordinates.
- Infinite leading model: RETURN090-B/091-C; reviewed Abel/variation
  identities RETURN100-D and true drift RETURN099-A; actual signed
  compensation and leading coefficients RETURN098-B/102-B.
- Common/augmented auxiliary Gram: RETURN092-B/C and RETURN088,
  formula attachment RETURN100-C; directed rational correction.
- Source assembly and terminal replay: RETURN091-B/C, RETURN092-A/C,
  RETURN093-B; even/odd fixed target 49/50 and 19/20.
- Error bands: RETURN095-A, RETURN098-B/099-A, RETURN102-A exact
  recomputation. The chart floor and common physical floor are
  also checked in RETURN099-B.
- Terminal algebra: RETURN100-B, including indefinite saddle trial
  residual and exact positive congruence transfer.
- Custody: RETURN100-A and RETURN102-C, 25 archives in immutable-by-
  policy distinct branches, original manifest 383 files.

## 5. Acceptance boundary and targeted remaining gate

The main remaining acceptance task is to bind the terminal
matrices from RETURN093-B to the actual regenerated source
assemblies RETURN091-C (even), RETURN092-A (odd) and reconstructed
auxiliary saddle RETURN092-C, by exact rational equality for TT/U/K,
rather than by the similarity of terminal margins. Also check that
the accepted independent retained Gram and literal moment records
are the corresponding parity/chart and anchors, and that the
latest independent error floors use the original margin and
comparator. This is an explicit cross-bank evidence gate, not
another search for favorable numerical values.

Historical component runs and analytic notes already support
the receiving argument, but there is not yet a single independently
replayed full dependency binding for all the supplemental evidence.
Record any mismatch and stop promotion. Successful component
receipts are not themselves hypotheses of RH. A self-contained
research-grade computer-assisted proof package still needs the
analytic derivations checked as a whole, the actual external
explicit formula citation checked if using the spectral statement,
and the documented numerical trust base. Native Lean formalization
is a separate goal, not logically required for a conventional
computer-assisted theorem. Do not conflate code execution with
a proof assistant's theorem certification.

No unconditional theorem is promoted in this note, and no
globalization argument has been supplied.
