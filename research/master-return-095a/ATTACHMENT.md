# RETURN095-A — Conditional original-form receiving lemma

This is an analytic interface, not a proof that the original Weil form
satisfies the interface. The numerical certificates apply to the frozen
diameter 5/2, retained real parity charts and specified comparator only.

## Constrained inverse identity

Let X be a real Hilbert space, A boundedly invertible self-adjoint positive
(or use its energy completion with bounded normalized row functionals).
Let C:X→R^c be the constraint rows, L:X→R^m the signed update rows.
Let J be symmetric invertible, not necessarily positive. Define
T=(C,L), J0=diag(0,J), K=J0+T A^{-1} T*.
Assume K invertible and the form D=A+L* J^{-1} L is coercive on ker C.
For t in the dual (identified via Riesz), put z=K^{-1}T A^{-1}t and
x=A^{-1}(t-T*z). Then T x=J0 z, hence Cx=0 and
Lx=J z_L. For every v in ker C,

    <D x,v> = <t-T*z,v> + <L*J^{-1}Lx,v> 
            = <t,v> - <z_L,Lv> + <J^{-1}Lx,Lv>
            = <t,v>.

Here the middle expression's update term is exactly
<L*J^{-1}Lx,v>; no positivity of K or J is required.
By uniqueness x is the constrained inverse response. Thus its energy is

    <t,x> = <t,A^{-1}t> - <T A^{-1}t,K^{-1}T A^{-1}t>.

For a retained source map t(p) this is TT - U* K^{-1} U with precisely
TT=t* A^{-1}t and U=T A^{-1}t in the same chart. This identity must be
attached to the actual constraint and signed updates, including the eight
restoration rows and the revised piecewise diagonal. In particular,
invertibility of the finite saddle alone does not imply positivity of D.

## Source perturbation and Schur transport

Suppose the original restricted quadratic form is bounded below by

    H[p] + D_phys[q] + 2 b(p,q),

where H is the positive retained form, D_phys is a positive complementary
energy, and its mixed-source dual norm is bounded by the modeled receiving
norm (with the stated 4/5 normalization). Suppose the model has squared
operator bound mu and the three source discrepancies, in the SAME
normalized receiving norm relative to H, have bounds e1,e2,e3.
Triangle inequality gives the full source squared bound
(sqrt(mu)+e1+e2+e3)^2. For r = sqrt(mu)+sum(ei)<1,
Cauchy and 2ab <= a^2+b^2 give

    original form >= (1-r)(H[p]+D_phys[q]).

This is conditional on all the displayed hypotheses, not a consequence of
finite matrix positivity alone. RETURN095-A only checks the rational
composition of the RETURN093-B modeled targets and RETURN085 error caps.

## Literal remaining obligations

1. Identify the actual form and retained/complement decomposition on a
   specified dense core; prove parity and constraint coverage, including
   excluded grid coordinates and pole-null condition if required.
2. Prove the comparator lower bound and its coercivity on the actual
   constrained complement, including signs and the 4/5 normalization.
3. Identify the raw retained matrix and the mixed source functionals with
   the original Gamma, prime and pole form; justify all infinite
   convergence and the cancellation/gauge identities.
4. Prove the equality of the modeled TT/U/K operators with the archived
   chart and reconstructed banks; attach the error reserves in the same
   inverse norm without switching comparator or domain.
5. Prove continuity/closed form extension for the original test class.
   A fixed diameter theorem is not an arbitrary support theorem or RH.

No native Lean proof, unconditional original-form positivity, arbitrary
support conclusion, or RH is asserted here.
