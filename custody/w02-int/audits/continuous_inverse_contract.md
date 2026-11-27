# W02 INT — Continuous inverse acceptance contract

Standing: analytical transport lemma checked; production result pending. No promotion of the 32-dimensional trial.

## Inputs and scope

- Frozen receiver: a4b63576f4a617c42975fdae85898fd1842739d6.
- Exact receiver SHA-256: 834d894be2bdac4bf7851d18862b203e7b234a68b6116f7d7d00a9c781929dac.
- Projection refinement: d06bdd70305b07195fce2fdf64e4945c0a130b62, run 46067568612.
- Independent arithmetic replay: 00107391ff270631dbaaf2b5b2eed469968506cd.
- Inherit BK/BS physical complement and signed operator hypotheses, including mu=53/10000 and M=I+B+R Delta^{-1}R* >= I on the physical Hilbert space. Do not infer circle-wide positivity of the signed symbol.

## Exact transport

Use columns V for the frozen physical receiver, U for any physical inverse response, E for the actual receiver discrepancy, and r=V-MU. All statements below are matrix Loewner inequalities, equivalently scalar inequalities for every column combination.

The variational identity is exact:
V* M^{-1} V = V*U + U*V - U*MU + r*M^{-1}r.

Since M>=I, r*M^{-1}r <= r*r and ||M^{-1}r|| <= ||r||.
For delta>0, t>0,
(V+E)*M^{-1}(V+E)
<= V*M^{-1}V + delta V*M^{-2}V + (1+1/delta) E*E.
This follows by Young on the cross term between M^{-1}V and E, and E*M^{-1}E <= E*E.
Furthermore M^{-1}V=U+M^{-1}r, so
V*M^{-2}V <= (1+t)U*U+(1+1/t)r*r.

Consequently if Fup >= V*U+U*V-U*MU, Gr >= r*r, Gu=U*U and D>=E*E, the complete upper comparison is
H = Fup + (1+delta*(1+1/t))*Gr + delta*(1+t)*Gu
plus (1+1/delta)*D.

There is no missing (1+delta) multiplier on Fup. This is not the cruder Young inequality in the M^{-1} metric.

For block errors E=sum E_i X_i, ||E_i||<=a_i, positive weights w_i summing to one give D=sum a_i^2 X_i*X_i/w_i by weighted Cauchy-Schwarz. Thus the existing penalty P=(1+1/delta)D/mu is sufficient. The original physical signed-energy lower matrix, minus P, H/mu and directed energy allowance, must be strictly positive. The exact head Gram defect <1 establishes dimension 32. No addition to earlier overlapping odd 31.

## Exact cell-projection identity

Let the doubled circle be [-4,12) and N equal cells. For real coefficient columns c_j, distinct positive labels j, f=sum c_j sin(j*pi*x/8)/2. Set d_k=sum_{j mod N=k} c_j exp(i(-j*pi/2+j*pi/N)) sinc(j*pi/N). The averages are (1/2) Im(sum_k d_k exp(2*pi*i*k*l/N)). Write d=A+iB. Discrete Parseval with cell measure 16/N gives
Gram(P_h f)=sum_k [(A_k-A_{-k})*(A_k-A_{-k}) + (B_k+B_{-k})*(B_k+B_{-k})].
Gram(f)=2 sum_j c_j*c_j. Their difference is the exact full-circle cell variance. Restriction to physical cells only decreases this positive Gram. This retains aliased and reflected labels, including j divisible by N, without spectral truncation.

## Executed local exact controls

SymPy symbolic direct cell integrals on N=4,8 with positive labels [1,N+1,2N+1,N-1,2N-1,N,2N,3N+2], two integer columns: exact zero discrepancy for the matrix identity. Additionally, a noncommuting rational metric M=I+B^T B for B=[[1,2,0],[0,1,1],[1,0,1]], V=[[1,2],[3,-1],[1,0]]/7, U=[[1,1],[0,-1],[2,1]]/11, E=[[1,-1],[2,0],[0,1]]/53, delta=3/100, t=7/8: both leading principal minors of upper bound minus exact inverse Gram strictly positive. These are local identity tests, not production execution or independent human review.

## Acceptance boundary

The production directed primitive/FFT/shift checks and exact positive terminal comparison are still required, followed by independent arithmetic replay of the immutable output. These do not newly prove inherited BK/BS operator identities or core transfer. A nonpositive lower comparison establishes insufficient bounds only. Until acceptance the conditional index stays <=205; prospective accepted 32 odd and 4 even yield <=204.
