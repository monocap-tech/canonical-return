# CR-INT-T1D-AO — Dyadic graph retention and chain-local Neumann bracketing

**Disposition:** analytical lower comparison, executed exact integer/Fraction inertia audit PASS. Conditional on the original fixed-window Gamma/Weil identity and AN's resolvent order; not Lean, remote CI, independent review, full positivity or RH. No earlier low-mode mixed certificate is subtracted from this bound.

## 1. Source and exact partition

Use I=(-2,2), the physical L2 metric, original doubly pole-null smooth core D^0, q0=Gamma-A. Split A=A2+Ar into the five prime powers 2,4,8,16,32 and the other 19 actual prime powers. All 48 signed translations are included. Put s=log 2, r=4-5s, so 0<r<s. In y=x+2, divide (0,r) into 48 equal phase intervals and (r,s) into 16 equal phase intervals. Translate each first-family phase interval by k*s for k=0,...,5 and each second-family phase interval by k*s for k=0,...,4. The 368 cells cover (0,4) up to endpoints; each phase chain consists of equal-length cells. Every admissible dyadic shift maps an entire cell to another cell of its chain (or outside the physical window), with unchanged local coordinate. In particular, no boundary fragment of a dyadic shift is dropped. A chain has length ell=r/48 or ell=(s-r)/16.

## 2. Exact source comparator

Let c_(j,n) be the vector of coefficients in the normalized Neumann cosine bases in the cells of chain j. Its squared Euclidean norm summed over j,n equals the physical L2 norm squared. The dyadic operator is exactly the direct sum over j,n of real symmetric matrices T_j with entries kappa_(2^|k-l|) for 1<=|k-l|<=5, zero diagonal. This follows by coordinate-preserving translation, with no approximation of log 2 by a rational grid.

For the residual operator Ar use AM's exact integer positive weight h (16384 cells). Weighted Schur gives <Ar f,f> <= integral (Ar h/h)|f|^2; equivalently this drops only the nonnegative residual weighted graph deficit. Use the 38 residual signed shifts, their certified AM shift floors and upper weight grids, and max of the two possible destination values (zero outside the window). Upward-round the resulting ratio on every AM cell to millionths. Take the maximum over every AM cell intersecting each new exact dyadic cell, using directed rational endpoint bounds to include a superset if necessary. Call the resulting diagonal potential v_j. This bounds all non-dyadic translations including those between chains.

Let C_j=diag(v_j)+T_j^+ +diag(e_j), where T_j^+ has upper grid weights U_k/10^9 and e_j is the sum of (U_k-L_k)/10^9 over incident retained edges. Exact rational log and integer sqrt enclosures verify L_k<=10^9*kappa_(2^k)<=U_k and U_k-L_k<=1. The row-sum error inequality via 2|ab|<=|a|^2+|b|^2 implies diag(v_j)+T_j <= C_j in quadratic-form order; merely replacing an off-diagonal entry by an upper bound would not suffice.

AN's variational direct-sum Neumann resolvent bound and complete digamma series imply

 q0[f] >= sum_(j,n) c_(j,n)^* [m_Gamma(n*pi/ell_j) I - C_j] c_(j,n).

The exact nonnegative terms omitted are the Gamma boundary matching defect, residual graph deficit and potential/retained-weight enclosure slack. The latter slack is positive semidefinite, not necessarily a scalar multiplication. This is a lower comparison, not an equality with the original form.

## 3. Executed exact inertia certificate

Use AN's improved lower digamma floor: x=n*pi/(2*ell_j), t_k=k+1/4, T=129/4. For x>T, m_Gamma >= m0+sum_(k<32)1/t_k - sum_(k<32)t_k/x^2 + log(x/T)+1/(2T)-T/(2x^2). This follows from convex trapezoid comparison of the decreasing positive digamma summand and its integral. The script uses a lower pi, upper ell, lower m0, rational lower log, and rounds the result down to billionths. n=0 uses the lower m0 directly. Every n>=1 used has x>T; the analytic lower expression increases with x. At margin eps=1/20000, compute exact LDL inertia of [gamma_lower(n)-eps] I-C_j, starting at n=0 and ending at the first positive definite matrix in each chain. Every pivot is nonzero; no float eigenvalue or floating-point inertia is used. Monotonicity of the full Gamma multiplier gives positivity for all higher n. The total negative inertia counts at the inspected modes bound the codimension of a subspace with q0>=eps||f||^2, by exact spectral decomposition of the finite symmetric comparator blocks. Restrict to the original smooth pole-null source and use the dimension intersection argument to bound its negative index.

Executed results: 64 chains, 368 cells; maximum first fully positive local mode 69; total exceptional dimension 7464, down 14110 from AN's 21574. The exact reflected chain matrices are checked to be related by row/column reversal, and phase chains occur in distinct reflected pairs. Reflection also multiplies the n-th cosine by (-1)^n. Thus each parity has exceptional dimension 3732. Conditional on the analytical imports, n_-(q0;D^0)<=7464. This is an exceptional basis of coupled chain modes, NOT a global Fourier cutoff at 7464. The previous AL mixed certificate cannot be subtracted. No actual negative direction was certified.

## 4. Next gate

AP should retain further prime shifts as exact operators on the phase variable (translations by log p modulo log 2 with appropriate chain-index changes). This must preserve physical metric and edge restrictions; the shifts for distinct prime bases are generally incommensurate, so a finite exactly aligned common cell partition cannot simply be assumed. An alternative is to bound the residual graph deficit on the finite exceptional chain space with rigorous mixed control. Restore Gamma matching defect in either route. A scalar norm estimate or positive compression alone cannot settle the missing spectrum.

Reproduce: `python3 ao_graph.py` with included AM dependencies (NumPy used for the supersolution proposal, all accepted inequalities independently checked with Python integers/Fractions). `audit.log` contains the executed output. Complete fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain OPEN.