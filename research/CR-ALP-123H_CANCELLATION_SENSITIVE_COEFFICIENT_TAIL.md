# CR-ALP-123H — Cancellation-sensitive coefficient tail envelope

Research-only; 123G joint trial bank 45921071905 still running at the start of this derivation. This derives a new infinite tail norm enclosure for the actual auxiliary residual of 123F, rather than recycling the absolute 5e-9 grouped-proxy error. No new numerical certificate or sign asserted.

## 1. Exact raw coefficients

Use t=x+a, a=11/10, v_j=j*pi/(2a), b=cosh or sinh in the corresponding parity, c_j=<e_j,b>=C_b v_j/(v_j^2+1/4), C_b=2b(0)/sqrt(a) in shifted coordinates. Write the normalized source f=sum_{i<=85} d_i e_i+F_b b. The joint trial q=-y=sum_{i in high bank} q_i e_i-tau h, h=(I-S85)b, tau=sum q_i c_i/||h||^2. Then f+q=sum t_i e_i+T_b b, with t_i=d_i+tau c_i for i<=85, t_i=q_i on the high bank, T_b=F_b-tau. Missing indices have zero coefficients. All equalities are symbolic, not claims that interval overlaps prove exact constraints.

Define A_i^G=v_i[R(v_i)-(-1)^i sum_m exp(-2a beta_m)/(beta_m^2+v_i^2)] and A_i^P=sum_n omega_n sin(v_i log n). These are the original 046 coefficients split before adding. Define L_i=d_i A_i^G+t_i A_i^P. For j greater than all trial indices, the actual raw high coefficient of Gf+P(f+q)+Dcap q is

    raw_j= -(2/a) sum_i L_i v_j/(v_j^2-v_i^2)
          +(2/a) A_j^G sum_i d_i v_i/(v_j^2-v_i^2)
          +(2/a) A_j^P sum_i t_i v_i/(v_j^2-v_i^2)
          +F_b G_bj+T_b P_bj-tau Dcap_j c_j,

where G_bj=Hb[j]-Pb[j] and Pb[j] is the literal signed prime-pole coefficient. This follows directly from 046 Hss and Hb, and from 123E's capped diagonal projection. On the pole-null Q85 tests the original pole term drops. The actual auxiliary error e=Q85(Gf+P(f+q)+Dcap q) has coefficient e_j=raw_j-alpha c_j for j>85, alpha=<h,Gf+P(f+q)+Dcap q>/||h||^2.

Alpha is computable from the same complete anchor data already used in 123E: <h,Gf>=sum d_i(hs[i]-ps[i])+F_b(hb-phb); <h,P(f+q)>=sum t_i ps[i]+T_b phb; <h,Dcap q>=sum q_i Dcap_i c_i-tau weighted. Here hs,hb are the full signed (Gamma+prime) pairings with h, ps,phb are prime-only, and weighted=<h,Dcap h> includes the complete infinite cap via Parseval. In particular alpha does not require a new integral or a sine tail truncation. Use the exact interval computations of 123E. Retain alpha separately; a finite projected sine row is not by itself a full norm.

## 2. Signed high-frequency cancellation

For j>J>=2048, Dcap_j=cap. Put A0=sum_i L_i and H0=-(2/a) A0-(tau cap+alpha) C_b. Combine these two leading nonoscillatory terms BEFORE taking an absolute value. Since c_j=C_b/v_j-C_b/(4 v_j(v_j^2+1/4)), the combined leading remainder is bounded by |(tau cap+alpha) C_b|/(4 v_j^3).

For any finite signed list x_i and integer p>=1, define M_l=sum x_i v_i^(2l), l=0,...,p-1, and H_p=sum |x_i| v_i^(2p). For v>=V>max v_i, the exact geometric identity gives

    sum x_i/(v^2-v_i^2)
       = sum_{l=0}^{p-1} M_l/v^(2l+2) + E_p(v),
    |E_p(v)| <= H_p/[v^(2p+2)(1-rho)], rho=(max v_i/V)^2<1.

Apply this separately to x=L, x=d_i v_i, and x=t_i v_i, retaining signed M_l as interval sums. In the L list the l=0 term has already been combined in H0. Use |A_j^G|<=4 and |A_j^P|<=4. The former follows from monotone lattice integral v sum_m 1/(beta_m^2+v^2) <=1/v+pi/4, and the reflected factor <=2; for v>2048 this is less than 4. The latter is the 024 sum-of-seven weights <4. No cancellation among unrelated unknown signs is presumed.

Safe pole column bounds for v>=1 are |G_bj| <= (80+2 log(1+v))/v and |P_bj|<=40/v. For the Gamma bound, use |g0|<6, |C_b|<3, |c_j-b0/(sqrt(a)v)|<=2/v, Lambda<=log(1+v)+9, and |Es|<=10/v^2 in 046's formula. The latter follows directly from sum |E_m|<10: |E0|<=2 sup|b|+2|b0|<6.4; for m>=1 use |E_m|<=(|b1|+|b0|/(4 beta_m)+exp(-2a beta_m)(beta_m|b0|+|b1|))/(beta_m^2-1/4); the first two sums bounded by (.6+.12)S1<.72, exponential term by the geometric tail (beta>=2.5, |b0|<1.2, |b1|<.6, denominator>=beta^2/2, exp(-a)<.4). Total <10. For the prime column use its exact expression, sqrt(a)>1, sum weights<4, |b(-a+log n)|<2, |b0|<1.2, |b1|<.6; the bound 40 follows. These crude pole bounds are multiplied by the actual signed, potentially very small F_b,T_b, never by the magnitude of the normalized source as a whole.

Combining gives a pointwise majorant for |e_j| as a positive sum of inverse powers of v_j and a log/v_j term. In particular the coefficient of v^-1 from the two large leading contributions is |H0|, not the sum of their separate magnitudes. The leading pole coefficients are not assumed to vanish; their saved intervals must decide whether this envelope is useful.

## 3. Full infinite norm with finite prefix

Use v_j>= (15/11) j (pi>3) and v_j< (10/7) j. For j>J, the sum of j^-s is <= J^(1-s)/(s-1), s>1, by decreasing integral comparison, even if all integers are counted instead of one parity. For the logarithmic component log(1+v_j)<=log(3j), and

    sum_{j>J} log^2(3j)/j^2
      <= [log^2(3J)+2 log(3J)+2]/J,

since this integrand decreases for j>=2048. Combine the resulting l2 upper bounds by Minkowski, using outward rational square roots and an upper log bound. This makes an executable rational norm enclosure for the full tail using only finite signed coefficient moments and elementary constants. Add the squared finite coefficient norm on 86..J (matching parity) to the squared infinite envelope. Divide by c=51/200 or use the 123F weighted bound; the cap-only part can instead use dcap=cap-21/8. All inequalities retain the original nonperiodic prime shifts.

This is a sufficient envelope, not proof that its margin passes. High-band finite coefficients up to J still need exact correlated evaluation. The original 046 psi_pair restricts imaginary/real shift ratio <1/3; extending the original generator unchanged to J=2048 would violate that guard. Use a larger shift or a validated alternative atan reduction, recheck Euler-Maclaurin/Taylor error guards, and reauthenticate the 85/165 prefix by interval overlap. This is an explicit implementation blocker. Do not claim bit-identical banks after changing the elementary evaluation algorithm.

## 4. Verification boundary

Classical analytical derivation only, under the same original source/interval and comparator premises. No actual H0, F_b, T_b or norm bound computed from 123G yet; no new runner, no sign or RH promotion. Next executable gate: compute alpha and signed moments from the preserved rational witnesses, verify the finite coefficient prefix and the high-frequency envelope separately, and then decide whether refinement is necessary.
