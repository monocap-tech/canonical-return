# RETURN108-K — Abstract compact-kernel and Vandermonde receivers staged

Status: five actual Lean theorem candidates staged in an isolated
Weil-Lab branch, not compiled or axiom-certified. No new workflow
activated. No change to the existing RECON staging branch or its
certified inventory. RH unresolved.

## 1. Source custody

Reopened RETURN108-J at aa0b727853f4a588bce90c28bfe4ee8948d91b3d.
RECON staging base: 486882068fa495d035e9d2450de0e4f215bf6e6c.
Toolchain and mathlib pin: v4.34.0 (literal lean-toolchain and
lakefile.toml inspected). The local container has no lean, lake or elan.

Created Weil-Lab branch research/return108k-abstract-kernel-staging
from that exact base, with new mathlib-only module
WeilDefect/MasterReturn/KernelCompleteness.lean.
Commit d83596db4c12ca9826b8b5168692d19a0ddaabc0,
blob 5661627b07544ca200aac74811292e8debb2dff3, reopened after creation.
Five theorem bodies and five exact #print axioms commands.
Static token scan found no sorry/admit/axiom admissions in the module.
This is not kernel verification. The Actions query for that head
returned no runs. No workflow file added or activated.

Mathlib v4.34.0 source inspected:
- Analysis/Normed/Operator/Compact/Basic.lean:
  IsCompactOperator.restrict, image_closedBall_subset_compact.
- Analysis/Normed/Module/FiniteDimension.lean:
  FiniteDimensional.of_isCompact_closedBall.
- LinearAlgebra/Vandermonde.lean:
  Matrix.eq_zero_of_forall_index_sum_pow_mul_eq_zero.
- LinearAlgebra/Dimension/Finite.lean:
  LinearIndependent.fintype_card_le_finrank.
RECON V08J/V09B give physical carrier lift/uniqueness and Mellin
entireness/uniqueness; source inspected but not imported by new module.

## 2. Actual staged declarations

Namespace WeilDefect.MasterReturn.KernelCompleteness:
1. finiteDimensional_of_compact_fixed:
   any closed submodule fixed pointwise by a compact operator
   is finite dimensional, by restriction and Riesz compact ball.
2. fullKernel_finite_of_compactSigned:
   for continuous linear T:H->E and R:E->H, R T=I+K with
   K compact implies ker T finite dimensional, since -K fixes it.
   No signed positivity or adjoint construction required.
3. independent_of_vandermonde_observations:
   vectors v_i with linear observations L_j(v_i)=c_j z_j^i,
   c_j nonzero and z_j distinct, are independent. Apply each row
   to a putative linear relation and cancel c_j; invoke the
   native Vandermonde zero-vector theorem.
4. zero_of_arbitrarily_large_kernel_banks:
   finite-dimensional kernel and arbitrary-size independent banks
   in that kernel conditional on x!=0 imply x=0.
5. fullSample_zero_of_enlarged_compact:
   composition of (2) and (4), abstract exact-kernel receiver.

The theorem names differ slightly from 108-J's initial suggestions:
we separated finite observation algebra from the analytic translate
attachment. The final receiver is a zero-vector implication with
explicit independent-bank hypothesis, not yet an actual sampling
injectivity theorem. The current source candidate may need Lean
elaboration/API repairs; no native compilation possible locally.

## 3. Simpler actual translation bank supplier

For nonzero supported f with entire Mellin F, choose z0 with
F(z0)!=0, using Mellin uniqueness. By continuity there exists
r>0 such that F(z)!=0 for |z-z0|<r. For any n>=1 choose d>0
with (n-1)d<b-a, and epsilon>0 with (n-1)epsilon<r.
Set s_i=i d and z_j=z0+j epsilon, i,j=0,...,n-1.
For the actual translated vector v_i=tau_(i d) f in H_b,

    F_(v_i)(z_j)=F(z_j) exp(d z_j)^i.

The nodes exp(d z_j)=exp(d z0) exp(d j epsilon) are distinct:
the first factor is nonzero, and the second is a strictly
increasing positive real sequence. Row factors F(z_j) are
nonzero. Thus the Vandermonde receiver supplies independence
using only continuity and the exact Mellin translation law,
without differentiating an exponential polynomial or invoking
complex identity theorem. For n=0 the empty bank is trivial.
Mellin evaluation at each z_j is a complex-linear functional on
the supported log carrier by finite-window integrability.
If the original full samples vanish, the literal translation
identity makes each bank vector lie in ker T_b; after lifting
the vectors to the kernel subtype, independence holds there too.
This is a paper-level supplier plan, not included as a compiled
Lean declaration in the staged module. General distinct-translate
independence is not needed for this contradiction.

## 4. Next

Native verification of the isolated five-declaration module;
then implement the actual translated physical lift, Mellin
linear observations and arbitrarily large kernel bank on one
fixed enlargement. The Green energy pairing, Dirichlet-to-log
density and graph-top closure transfer remain separate.
No arithmetic premises removed from the signed compact-defect
supplier; positive-channel coercivity is still a different
obligation. No proof of RH or mathematical positivity.
