# CR-ALP-062 candidate — four-moment matrix inverse comparison

DRAFT / NONPROMOTING. Parallel continuation of CR-ALP-060, not a replacement for the CR-ALP-061 physical integration. No new physical moments are claimed.

Let R map the finite source space into Q85, L=gamma I+P_Q with 0<ell I<=L<=u I, gamma=(ell+u)/2, s=(u-ell)/2. Define G=R*R, M=R*P_Q R, N=R*P_Q²R, T=R*P_Q³R and J=R*L^-1 R. All moment matrices are Hermitian. Write H=(L-ell I)(u I-L), H0=R*HR=s²G-N and H1=R*LHR=gamma H0+s²M-T.

For S=L^-1/2 H^1/2 R and W=L^1/2 H^1/2 R, S*S=gamma G-M-ell u J, S*W=H0 and W*W=H1. For every source matrix Y, expanding (S-WY)*(S-WY)>=0 gives

J <= (gamma G-M-H0Y-Y*H0+Y*H1Y)/(ell u).

Therefore ell u K-gamma G+M+H0Y+Y*H0-Y*H1Y>0 suffices for full Schur positivity. Optimize by Y=H1^dagger H0, since ker H1 is contained in ker H0 by the Gram factorization. No inverse is required to check a rational frozen trial. Setting Y=a^-1 I recovers the scalar cubic majorant, and taking one source column gives the CR-ALP-061 sharp scalar upper bound. No claim of general matrix moment extremality beyond this trial optimization.

For the lower comparison, use X=RA+P_Q RB for arbitrary source matrices A,B. Then J>=R*X+X*R-X*LX by the positive L-norm square of L^-1R-X. With D=[G;M], Z=[A;B] and block F=[[gamma G+M,gamma M+N],[gamma M+N,gamma N+T]], the lower bound is D*Z+Z*D-Z*FZ. Optimize by Z=F^dagger D, since F is the L-Gram matrix of [R,P_Q R] and ker F is contained in ker D*. Again arbitrary rational frozen trials require no inverse certificate. Lower comparison exceeding K in a direction rejects this sufficient L comparison, not the original Weil form.

Exact BigInt rational regression executed in the tool runtime: arbitrary noncommuting rational trials on a four-dimensional spectral example, exact Gaussian lower equality, exact Lobatto upper equality, and singular endpoint equality all PASS. Code: src/audits/cr_alp_062_matrix_four_moments.js (commit 91b29a95706a0b3b860129d910b845cc13d985f0). Replay with node. This regression is finite abstract arithmetic, not native Lean or physical integration.

The full coupled matrix moment banks are not yet evaluated. Two scalar witnesses cannot establish a full matrix inequality. If scalar CR-ALP-061 fails to decide, this comparison offers an alternative to blindly raising the scalar moment order, but computing and transporting the full banks remains a separate obligation. Keep this on the CR-ALP-060 candidate branch while 061 runs. Inputs and new source remain repository-controlled; historical standalone archive reconciliation remains open. No receiving-system promotion.
