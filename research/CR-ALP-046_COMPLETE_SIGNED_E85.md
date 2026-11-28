# CR-ALP-046 — Complete original signed E85 finite section

Status: exact local arithmetic PASS on the finite section, conditional on the classical source identities. No Lean compilation, mixed E85/Q85 Schur certificate, whole-collar positivity, arbitrary-support transfer, or RH claim. The earlier CR-ALP-045 Q85 coercivity bound 1/4 is unchanged.

## Source and metric

On I=(-11/10,11/10), let e_j=a^(-1/2) sin(nu_j(x+a)), nu_j=j*pi/(2a), a=11/10. Let b=cosh(x/2) in the odd-index sine sector and b=sinh(x/2) in the even-index sector. c_j=<e_j,b>, B=<b,b>. Project by v_j=e_j-c_j b/B. The actual physical Gram is M_ij=delta_ij-c_i c_j/B in each sector. The equivalent predecessor basis consists of w_j=v_j-(c_j/c_k)v_k for j>=3, and u_k=v_k for k=1,2. Both reflection sectors are orthogonal for the signed form. The matrix is restricted to both original pole constraints; it is not an unrestricted form.

The attached source is the same CR-ALP-029/043 analytic expression for the signed prime-Gamma form: g(r)=Re psi(1/4+ir/2)-log pi and the genuine negative paired shifts at n=2,3,4,5,7,8,9 with omega_n=Lambda(n)/sqrt(n). No approximation of the original signed prime set or removal of the Gamma series.

## Full Gamma entries without long resolvent truncation

Put L=2a, beta_m=2m+1/2, sgn=(-1)^j, s=-sgn, and R(v)=sum_m 1/(beta_m^2+v^2). CR-ALP-043 established R(v)=Im psi(1/4+iv/2)/(2v) and the off-diagonal sine entry. Differentiating this identity at v>0 gives

    R2(v)=sum_m 1/(beta_m^2+v^2)^2
         = Im psi(1/4+iv/2)/(4v^3)-Re psi'(1/4+iv/2)/(8v^2).

The differentiated series converges locally uniformly. With R_s(v)=R(v)-sgn sum_m exp(-L beta_m)/(beta_m^2+v^2), and Theta(v)=sum_n omega_n sin(v log n), set A_i=nu_i R_s(nu_i)+Theta(nu_i). Then the complete signed sine matrix in matching parity is

    H_ij=-2/a*(nu_j A_i-nu_i A_j)/(nu_j^2-nu_i^2), i != j;
    H_ii=g(nu_i)-2 nu_i^2/a*[R2(nu_i)-sgn sum_m exp(-L beta_m)/(beta_m^2+nu_i^2)^2]
           -sum_n omega_n*[(L-log n)/a*cos(nu_i log n)+sin(nu_i log n)/(a nu_i)].

The diagonal identity follows from the CR-ALP-029 exact resolvent diagonal and the sum 2 sum nu^2/[beta(beta^2+nu^2)]=Re psi(1/4+i nu/2)-psi(1/4). Off-diagonal identity is the uncentered CR-ALP-043 formula. All seven prime terms are exact analytic expressions.

## Both genuine boundary-anchor columns

Use t=x+a, b(t-a)=b0 cosh(t/2)+2b1 sinh(t/2), b0=b(-a), b1=b'(-a). For v=nu_j, c_j=2 b0 v/[sqrt(a)(v^2+1/4)], Lambda(v)=Re psi(1/4+iv/2)-psi(1/4), and E_m=int_0^L exp(-beta_m t)b(t-a)dt-b0/beta_m. The uncentered Gamma sine-pole entry is

    g(0)c_j + Lambda(v)[c_j-b0/(sqrt(a)v)]
              -2v/sqrt(a) sum_m E_m/(beta_m^2+v^2).

The E_m divided differences are exactly those in CR-ALP-043 section 2. The paired prime translation contributes -omega_n*2/[sqrt(a)(v^2+1/4)]*[v(b(-a+log n)+b0 cos(v log n))-b1 sin(v log n)]. No pole-aligned piece is dropped before projection.

For the missing b-b Gamma entry, use K_beta with kernel exp(-beta|x-y|), and J_beta=2I/beta-K_beta. For beta>1/2 solving the elementary resolvent differential equation with K'(0)=beta K(0), K'(L)=-beta K(L) yields

    K_beta b(t)=2 beta b(t)/(beta^2-1/4)
       -(beta b0-b1)/(beta^2-1/4) [exp(-beta t)+s exp(-beta(L-t))].

Also int_0^L exp(-beta t)b(t)dt=[beta b0+b1-s exp(-L beta)(beta b0-b1)]/(beta^2-1/4). Hence for m>=1 the J_beta quadratic value is

    -B/[2 beta(beta^2-1/4)]
    +2*(beta^2 b0^2-b1^2)/(beta^2-1/4)^2
    -2*s exp(-L beta)*(beta b0-b1)^2/(beta^2-1/4)^2.

Set S1=sum_{m>=1}1/(beta_m^2-1/4)=1-log 2, S2=sum_{m>=1}1/[beta_m(beta_m^2-1/4)]=2psi(5/4)-psi(1)-psi(3/2), S4=sum_{m>=1}1/(beta_m^2-1/4)^2=pi^2/6-1-2S1. The summed b-b Gamma energy is

    g(0)B + J_(1/2)(b,b) - B S2/2 + 2 b0^2 S1
       +2(b0^2/4-b1^2)S4 - 2s sum_{m>=1} exp(-L beta_m)*(beta_m b0-b1)^2/(beta_m^2-1/4)^2.

At beta=1/2 use the exact triangular double integral of the two exponentials exp(+-t/2), with the resonant term int t exp(-t)dt or int t dt. The b-b prime entry is -2 sum_n omega_n int_0^(L-log n) b(t)b(t+log n)dt, evaluated by the same exponential primitives. Finally project the full matrix by v_j=e_j-c_j b/B and change to w_j,u_k. The entire E85 block has dimensions 43 and 42.

## Interval implementation and error controls

Use integer intervals at scale 10^340 with outward rounding. Machin pi, atanh log, arctan, Taylor sin/cos and exponential use 400-term series with rational error comparisons. For psi and Re psi', shift the argument by 512 and use 90 Bernoulli terms. The classical integral remainder for psi, using sup|B_180({t})|<=4*180!/6^180, has absolute bound 4*180!/(6^180*180*512^180); differentiating under the integral gives bound 4*180!/(6^180*512^181) for psi'. Both are <10^-290, checked exactly, and each returned interval is widened by that radius. For the shifted arguments the arctangent ratio is <1/3. The initial 120-decimal attempt suffered amplification by large Bernoulli coefficients; it provided no certificate. At 340 decimal places the row errors are <10^-150 before outward compression. All exponential corrections use 256 terms, with omitted tail bounded by 100*(2/5)^(4*256+1)/(1-(2/5)^4)<10^-340. The factor 100 dominates the rational coefficients of the omitted terms: beta>512, nu>=1, |b0|<6/5, |b1|<3/5, beta^2-1/4>beta^2/2. The coefficient denominators are positive; all elementary quantities are enclosed by the interval routines.

The builder checks 85 physical pivots after subtracting 10^-110 M and an outward row-radius bound from each integer midpoint matrix. An independent verifier compresses all entries outwards to 70 decimal places, checks positive definiteness of both physical Gram blocks, and performs a separate lower-triangle fraction-free determinant algorithm. With 70-digit intervals the row errors are <10^-65. The exact shifted tests establish 43/43 pivots with shift 10^-40 (odd indices) and 42/42 pivots with shift 10^-30 (even indices). Larger trial shifts were inconclusive, not certified negative original witnesses. The common finite margin is 10^-40. Exact predecessor regression checks 98 physical Gram overlaps, positive diagonal increments over the old E14 lower matrix, and entrywise compatibility with the Cauchy-Schwarz bounds of its omitted positive Gamma remainder. That is not independent reconstruction of the new entries or a general PSD test of the difference.

Matrix SHA-256: d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22.

## Result and next obligation

Conditional on these continuum identities, W(p,p)>10^-40 ||p||^2 for nonzero p in E85. CR-ALP-045 separately proved W(q,q)>1/4 ||q||^2 for nonzero q in Q85 H0 in the logarithmic form domain. The two subspaces together span H0, but their mixed signed interactions have not been bounded: separate positivity is insufficient. The exact sufficient scalar Schur gate is H85-4 B85 PSD with B85 the full projected residual Gram, or a stronger frequency-dependent inverse comparison on the same source. Neither has been evaluated here. The original nonpositive-index bound 66 is not improved just by the finite result. No formal or receiving acceptance ledger promotion.

Reproduction: run python cr_alp_046_full_finite.py (substantial integer computation), then python check_saved.py and python check_predecessor.py. The independent verifier can run directly on the included pinned bank. It checks the interval matrices rather than rebuilding analytic inputs.
