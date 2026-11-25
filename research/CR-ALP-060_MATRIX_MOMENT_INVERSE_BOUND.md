# CR-ALP-060 candidate — matrix-valued three-moment inverse bounds

DRAFT / NONPROMOTING. Parallel candidate prepared while CR-ALP-059 native integration runs. It does not require or assert a new numerical second moment, full matrix bank, original Weil positivity or Lean result. Parent pin: 4b60bff3e0bde60d186b216741f5a2f68e15c718. The separate branch avoids changing the running 059 verification commit.

Let R map the finite source space into Q85, L be bounded selfadjoint with 0<ell I<=L<=u I, J=R*L^-1 R, G=R*R, M1=R*LR, M2=R*L^2R. For W=L-ell I let S=L^-1/2 W^1/2 R, T=L^1/2 W^1/2 R. Then S*S=G-ell J, S*T=B=M1-ell G, T*T=A=M2-ell M1. Expansion of (S-TY)*(S-TY)>=0 gives, for every finite matrix Y,

J <= (G-BY-Y*B+Y*AY)/ell.

In particular, ell K-G+BY+Y*B-Y*AY>0 suffices for K-J>0. Taking Y=a^-1 I recovers exactly the CR-ALP-058 scalar quadratic family. Since A>=0 and ker A is contained in ker B (T*x is not needed: Ax=0 implies Tx=0, hence Bx=0), the optimal trial correction is B A^dagger B, attained by Y=A^dagger B. A singular A is harmless in the arbitrary-trial form.

Likewise W=u I-L gives S*S=u J-G, S*T=D=u G-M1, T*T=C=u M1-M2; hence J>=(G+DZ+Z*D-Z*CZ)/u for arbitrary Z, optimized at Z=C^dagger D. The bounds are Loewner inequalities, and the optimal matrix upper bound is no weaker in each direction than any common scalar quadratic bound; no general claim of sharpness among all operator-valued moment measures is made.

In prime coordinates M1=gamma G+M and M2=gamma^2 G+2 gamma M+N where M=R*P_Q R, N=(P_Q R)*(P_Q R). All entries are linear combinations of the same three full moment banks. Rational Y,Z can be selected numerically and frozen, then the final congruences and shifted LDL verified with outward interval matrices. No certified inverse of a numerical A is necessary.

Exact abstract audit: full two-dimensional nonorthogonal source, L=diag(1,4), R=[[1,1],[0,1]], K=J+G/100. The matrix trial attains J exactly. A scalar interpolation a>0 has excess (t-ell)(1-t/a)^2/(ell t) on the two eigen-directions. Passing requires |1-t/a|<1/10 for both t=1,4, which is impossible since the necessary intervals for 1/a are disjoint. Proper residual range and arbitrary trials tested on a 3x2 example; singular A,C tested at ell,u. Standard-library Fraction audit locally PASS (three output checks), fetched repository source equals submitted text. This is an abstract regression, not an actual Weil witness.

Prerequisites to the full matrix gate: 059 actual two witness moments and their endpoint verification; if useful, evaluate all global G,M,N entries on the frozen source columns, then transport the proxy errors as a matrix norm or individual intervals. Positivity on two witnesses is not PSD on all 85 directions. Rejection of the L comparison is not a negative original Weil direction.

Custody: source audit and this note in Git, no raw data created or left only locally. CR-ALP-058 input archive and its existing custody inventory retained. Do not promote or merge this candidate into receiving-system ledgers before its own evidence checks. The historical standalone archive reconciliation remains open.
