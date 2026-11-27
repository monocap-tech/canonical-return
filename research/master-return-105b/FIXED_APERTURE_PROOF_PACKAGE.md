# MASTER-RETURN-105-B — Fixed-aperture computer-assisted proof package

Status: integrated candidate for mathematical acceptance, not a declaration
of unconditional theorem certification. The scope is a single aperture.
The final custody workflow RETURN105-A (45918163435) has passed its
literal scalar-attachment step; its archive-preservation step is separate.

## 1. Statement

Let a=5/4, I=(-a,a), and f be a real L2 function supported in the
closure of I. Put F(z)=integral_I f(x)e^(izx) dx and
c(t)=integral_R f(x)f(x+t) dx, using zero extension. Define

    W(f)=1/(2 pi) integral_R |F(r)|² g(r) dr
           -2 sum_{n in N} (log p)/sqrt(n) c(log n)
           +2 F(i/2) F(-i/2),

where g(r)=Re psi(1/4+ir/2)-log(pi), N={2,3,4,5,7,8,9,11},
and n is a power of the prime p. The domain is the supported
logarithmic form domain: integral log(2+|r|)|F(r)|² dr finite.

The asserted computer-assisted inequality, conditional on the
correctness of the listed analytic and arithmetic certificate premises,
is

    W(f) >= (3/10^59) ||f||².

Translation extends it to real tests supported on any interval of length
5/2; no dilation or arbitrary-support extension is asserted. The
Guinand–Weil formula identifies W with the sum of F(z_rho)F(-z_rho)
on smooth compactly supported real tests, where z_rho=(rho-1/2)/i
and zeros are counted with multiplicity. This identity is unconditional
and is a corollary, not a premise for proving the physical inequality.
The spectral summands are not assumed individually nonnegative.

## 2. Analytic decomposition

All three terms of W are invariant under reflection, hence even and
odd tests decouple by polarization. Work in one parity. The complete
interval grid uses r_k=2 pi k/5 for positive even k; include the constant
in even parity. Nonconstant modes have squared L2 norm 5/4. The shifted
odd-frequency grid is independently complete in each parity (modulation
of the complex basis); its nonconstant real modes have the same norm.

The finite retained space P has 87 physical generators in each parity:
81 positive even frequencies, five unpaid positive odd frequencies,
and the constant (even) or sinh(x/2) (odd). These labels and their
rational coefficient chart C/T0 are the authenticated RETURN059 and
RETURN096-A objects. The chart is nonsingular. Orthogonally decompose
f=p+q, with p in P and q perpendicular to P. Write p in chart
coordinates x and q in rephased complementary grid coefficients y,
so ||q||²=(5/4)||y||². The complementary grid J consists of the
16302 paid finite even coordinates through 32766 and all even k>=32768.
The five odd-frequency orthogonality equations (plus sinh in odd parity)
are bounded linear functionals of y. Their Newton divided-difference
row basis and invertible Z change preserve the kernel, denoted ker C0.
See RETURN096-A and RETURN100-C for the actual index/rank checks and
partial fraction identity.

The series g(r)=g(0)+sum_{m>=0} 2r²/[beta_m(beta_m²+r²)],
beta_m=2m+1/2, gives the line exponential-kernel form via Fourier
inversion. The nonnegative partial Gamma increments converge
monotonically. Periodize on a circle of length 5. The image difference
at x,y in I is
2 exp(-5 beta) cosh(beta(x-y))/(1-exp(-5 beta)). The sum of
image differences defines a bounded operator, by summing the geometric
norm majorant for cosh/sinh. This also proves equality of the line
and periodic logarithmic form domains. Prime translations do not
wrap, since log n<5/2 and 5-log n>5/2. Pole terms are bounded.
The retained generators have O(1/k) circle Fourier coefficients, so
belong to the form domain; thus both p and q belong. For fixed p the
log-multiplied Fourier coefficients remain in l2, so the mixed form
extends in q by L2 continuity. Partial Fourier sums need not satisfy
the constraints: take the limit first and impose C0 y=0 on the full y.

## 3. Complementary lower comparator

For even q the image and pole corrections are nonnegative. For odd q
the sinh constraint cancels the pole and the first negative image;
the m>=1 image cost is less than .001||q||² by the bound

    exp(-25/4)/[5(1-exp(-25/2))(1-exp(-5))] < 1/1000.

The original periodic Fourier symbol is

    s(k)=g(r_k)-2 sum_{n in N} (log p)/sqrt(n) cos(r_k log n).

Circle coefficients of zero-extended functions are half the interval
coefficients. Its form is (5/8) times the sum of the two interval-grid
symbol forms. The exact signed symbol check RETURN096-B pays .1
baseline in each grid plus its individual assigned surplus. The
finite signed ledger uses v_j-d_j+u_j, with no double counting of
restoration weights. The resulting lower comparator has diagonal
A_k=.099+w_k/2 on the finite paid grid and A_k=.174 on k>=32768.
In odd parity strengthen it to .519 on 65536<k<=262144 and
1.199 on k>262144. Those inequalities were checked against the
monotone lower symbol and a prime triangle bound.

For the Cauchy rows h_j(k)=j/(k²-j²) in even parity and
k/(k²-j²) in odd parity, the original odd-grid coefficient is,
up to sign, (4/pi)h_j(y). The signed energy inside the comparator
is therefore 8(v_j-d_j)|h_j(y)|²/pi². The rational row scaling s_j
cancels with its metric J_j=+/-s_j² pi²/(8 weight).
For the odd restoration bank B (8 by 500), R is positive and
R=B diag(pi_upper²/(8u_j)) B^T. The projection inequality
(Bh)^T R^-1(Bh)<=sum 8u_j h_j²/pi² justifies its lower comparator
contribution. Collect all signed and restoration rows into L,
with invertible block metric J_s. In the actual constraint kernel,

    W(q)>=(5/4)D[y],  D=A+L* J_s^-1 L.

The unprojected negative update G_- satisfies
||G_- A^-1 G_-*||<=13/25 by RETURN095-B, for both parities.
The bound persists under the odd strengthening. Dropping positive
updates gives D>=(12/25)A>=(12/25)(99/1000)I.
Consequently D_phys=(5/4)D has physical mass lower bound
(297/6250)||q||². No assertion that the auxiliary saddle itself
is positive is made.

## 4. Retained and mixed physical forms

H[x]=W(p) is the actual retained physical form in the rational chart.
The independently enclosed raw 87 by 87 matrices of RETURN101-C
are contained in the accepted RETURN093-B terminal matrices. The
separate retained scalar margin h and exact chart mass bound M satisfy

    H[x]>=h||x||²,  ||p||²<=M||x||²,
    M=(5/2)*87*||C/T0||_F².

These bounds are checked in RETURN099-B/103-B. They do not assume
the chart orthonormal.

For the mixed form write W(p,q)=<t(x),y>. For a single retained
generator p_j write ell=p_j''/p_j (constant eigenvalue),
b=p_j(a), v=p_j'(a), M_j=integral p_j exp(x/2), and r=2pi k/5.
The boundary differential identity gives the cosine/sine mass
2 sigma v/(r²+ell) or -2 sigma r b/(r²+ell), sigma=(-1)^(k/2).
RETURN098-C derives the full Gamma, prime and pole grouped rows,
including the beta=1/2 sinh resonance M_j=sinh(a)-a. The actual
raw coefficients (ell,b,v,c,d) are independently evaluated in
RETURN101-B, and their charted mixed entries in RETURN104-A.
The pole d uses its separate denominator r²+1/4. For the exact
resolvent vector moments B(x),V(x),L(x) of RETURN099-A, the
grouped mixed rows take the form

    t_even=alpha B/r+d_scalar V/r²+2 L/r²,
    t_odd=-d_scalar B/r+alpha V/r²+2 L/r.

Here alpha=2r(h_gamma-e_gamma)+2Ps and
d_scalar=g(r)+g(0)+4-1/(r²+1/4)+2f_gamma-2Pc.
h_gamma,e_gamma,f_gamma denote the convergent m>=1 resolvent
series, Ps,Pc the eight weighted prime sine/cosine sums.
No raw coefficient is confused with d_scalar.

Use four disjoint bands: 8111 explicit rows through 16384;
8192 model rows 16384<k<=32768; 114688 high-band model rows
32768<k<=262144; and every even k>262144. The odd reweighting
of 98304 high-band rows 65536<k<=262144 is a correction,
not an additional support partition. The finite models truncate
the geometric resolvent at N=12 and N=8 respectively, with
alpha0=pi/2+2Ps and d0=D0+log(k/K)-2Pc,
D0=g(2pi K/5)+g(0)+4. The odd infinite model includes the
actual compensated tail -G r/(r²+1/4), with G=D0 B0-2L0.
The odd leading scalar is (2Pc-g(r)+g(R))/r, not a pure log
model. The compensation is retained in TT, U and restoration
crosses; global sinh orthogonality does not cancel a tail-only row.

For the error bounds use the mesh integral variation estimate
|g'(r)-1/r|<3/r², hence the weaker inherited drift <=8/R,
and resolvent estimates |alpha-alpha0|<=6/r,
|d_scalar-d0|<=8/R+3/r². Uniform absolutely convergent geometric
series in (K/k)² and positive coefficient majorants give the
omitted moment bounds; the infinite odd compensation difference
has the explicit 1/(4r³) reserve. RETURN099-A gives the full
six envelopes. Their square sums and restricted inverse trace
are independently recomputed in RETURN102-A. The finite moments
are bound literally in RETURN102-B.

## 5. Constrained inverse and terminal payment

Put T_aux=(C0,L), J0=diag(0,J_s),
K=J0+T_aux A^-1 T_aux*, TT=m* A^-1 m, U=T_aux A^-1 m,
where m is the four-piece modeled mixed source. All auxiliary
rows are bounded l2 functionals. A and D are bounded and
coercive. The exact terminal inverse residual establishes
invertibility of the finite symmetric K.

For z=K^-1 Ux and y0=A^-1(m(x)-T_aux* z), one has
T_aux y0=J0 z. Thus C0 y0=0 and Ly0=J_s z_L.
Testing against any y in ker C0 gives
<D y0,y>=<m(x),y>. The squared physical dual norm is therefore

    (4/5) x^T (TT-U^T K^-1 U) x.

The factor 4/5 is the physical inverse of 5/4, not a chart
normalization. The exact trial identity with E=U-KY is
TT-U^T K^-1 U=TT-Y^T U-U^T Y+Y^T K Y-E^T K^-1 E.
The last term is bounded above by ||K^-1|| ||E||_F² I;
no sign of K^-1 is assumed. The terminal checker obtains an
inverse cap from the Neumann residual of the symmetric rational
trial inverse, and directed entrywise reserves for the trial
core. It uses a weighted raw diagonal radius majorant before
the rational H congruence. Nonsingular modular-rank-verified
congruences and either strict row dominance or exact Bareiss
positive pivots certify the comparisons. RETURN093-B obtains

    (4/5)(TT-U^T K^-1 U) <= mu H,
    mu=49/50 even, 19/20 odd.

The three error bands have physical dual norms bounded by
eps_i sqrt(H[x]) by the negative-update comparison and
restricted Woodbury trace. The full source norm is bounded by
r sqrt(H[x]), r=99/100+sum eps_i even or 39/40+sum eps_i odd.
Both rational square-root caps dominate sqrt(mu). The triangle
inequality, not orthogonality of error bands, is used.

## 6. Conclusion of the receiving implication

    W(f)>=H[x]+D_phys[y]+2<t(x),y>
        >=H[x]+D_phys[y]-2r sqrt(H[x]D_phys[y])
        >=(1-r)(H[x]+D_phys[y])
        >=(1-r) min(h/M,297/6250) ||f||².

The last mass equality follows from L2 orthogonality of p,q.
RETURN102-A and RETURN103-B check the rational floor in both
parities; decimal diagnostics are 3.3925169506380073e-59 even,
9.593091391347799e-55 odd. In particular both exceed 3/10^59.
Sum the orthogonal parity inequalities. The receiving implication
is complete conditional on the literal analytic/interval premises.

## 7. Reproduction and acceptance record

The original reproduction archives are reconstructed with
archive/master-return/reconstruct.sh; this verifies part and whole
SHA-256 checks. The 383-file source manifest is pinned at
38b4781ddac2b4f0437f75f8eb802e78a54c1fbc366f4690880aa73b0aa9d112.
Do not edit original archives. Later evidence archives are
preserved under distinct custody commits.

- RETURN096-A/B: index partition, actual constraint kernel,
  signed symbol, physical periodization and normalization.
- RETURN094-A/095-B: restoration metric and negative-update cap.
- RETURN098-A/B/C, 099-A, 100-B/C/D, 103-C: analytic derivations
  and receiving review, not proof-assistant declarations.
- RETURN101-A/B/C, 102-B, 104-A: Arb source enclosures, including
  2,001,138 finite source comparisons. RETURN105-A binds the
  actual high-band D0/pi endpoints and 900 direct moment entries.
- RETURN087 -> 104-C: finite and feature contraction banks
  reproduced by exact fmpz_mat A^T B. The actual even 091-C
  and odd 091-B assembly source explicitly pins these same
  RETURN087 banks and consumes their intervals. Neither reverts
  to the historical broader parent values. Exact bank equality
  removes the binary64 dot-product assumption from these inputs.
- RETURN088/089/090-A: exact rational/object-integer Gram and
  high-band contractions; scalar input enclosures bound above.
- RETURN090-B/091-C -> 104-D: Arb phase/chord/Gamma replay
  reproduces the actual consumed infinite-tail banks exactly.
- RETURN092-B/C -> 104-B: common/augmented auxiliary Gram
  exact fmpz_mat equality (815414 entries), four terminal inputs
  identical; baseline and strict-tail reciprocal corrections
  preserved.
- RETURN091-C (even), RETURN091-B -> 092-A (odd), 092-C
  (auxiliary replacement), 093-B (retained replacement):
  exact terminal lineage confirmed by 103-B. The older overlap
  checks are not used as a proof of source validity.
- RETURN102-A/103-B: six exact error budgets and two final
  physical floor rational checks.

The mathematical trust base is Python integers/Fraction, exact
NumPy object-integer operations, FLINT exact integer matrices,
Arb rigorous balls, and the analytic identities and inequalities.
Recent replacements pin python-flint 0.8.0 and numpy 2.2.6.
Reference readers may still import mpmath/SciPy; no claim of a
dependency-free runtime or zero trusted software is made.
An authenticated receipt certifies an execution, not the truth
of the analytic assumptions. All code must be executed without
disabling assertions. A proof assistant is not used here.

## 8. Terminal boundary and next action

Finish RETURN105-A custody and verify the remote commit and
manifest. Record the pinned exact backend and terminal lineage
receipt digests in the final reproduction index. Subject the
integrated theorem to independent mathematical review, with
particular attention to the completeness of the analytic
attachments, the selected interval arithmetic implementation
and the domain definitions. Do not label this internal assembly
as external referee acceptance. No additional numerical search
is warranted without an identified flaw or failed binding.

No unconditional theorem promotion is made here. The certificate
is a fixed-aperture physical coercivity implication with extensive
verified numerical premises. A global Weil criterion requires
arbitrary compactly supported tests, and no transport beyond
diameter 5/2 has been established. RH remains unresolved.
