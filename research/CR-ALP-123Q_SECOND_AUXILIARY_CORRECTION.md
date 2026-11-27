# CR-ALP-123Q — Second auxiliary correction and uniform polynomial inverse proposal

Research-only. The 123O expanded residual receiver (run 45925861807) remains executing. The 123P workflow is staged with explicit activation. This note develops a conditional alternative if its generalized matrix-order sufficient bound remains inconclusive.

## 1. Retain the coupling before applying worst-case inverse error

Use the 123F definitions with fixed normalized source f in the span of E85 and b, and y in Q85:
e=r-Cy, b_y=a-X*r-r_C* y, E_joint=W(f,f)-[2 Re<r,y>-<y,Cy>].
Take any second auxiliary trial w in Q85. Set y'=y+w. Then the exact identities are

    e'=e-Cw
    b_(y')=b_y-r_C* w
    E_joint'=E_joint-2 Re<e,w>+<w,Cw>.

If U' bounds <e',C^-1 e'>, the 123K/123P matrix order receiver yields

    sigma >= E_joint' - U'
      - ( ||S0^-1/2 b_(y')|| + sqrt(U'/t) )^2,

where S0^-1=K^-1+N^-1, N>=tK. Proof: substitute the exact identities into 123F and use the previously derived Loewner comparison. No need to charge the original e by the factor 1/t. A candidate w can be chosen in a much larger Q85 projected bank independently of the expensive ill-conditioned full joint low-block optimization. This correction need not help: the new finite coupling b_(y') may increase, and the energy may decrease. Both must be explicitly evaluated. No claim that this is equivalent to reoptimizing the low trial.

For w=sum w_j Q85 e_j, <e,w> uses the actual projected sine coefficients. r_C* w=R* w-X* Cw uses the original mixed Gamma/prime columns and the frozen correction bank. In particular Cw=e-e' gives X* Cw=X*e-X*e'. All these vector terms are signed and must be evaluated before taking norms. The full pole bridge is retained in w; high sine projection alone would change the operator. The required new data are the rational w coefficients, signed Cw and mixed finite pairings, and complete norm upper bound on e'. 123O currently saves X*e but does not save its entire explicit residual prefix, so a new coefficient evaluation or source regeneration is needed to construct w.

## 2. A bounded positive operator inverse proposal, independent of the near-null low block

The actual common capped C=Q85(Dcap+P)Q85 satisfies

    ell I <= C <= u I, ell=51/200, u=cap+21/8=7.072938456.

The lower bound uses 068/122D and the exact prime norm bound; the upper follows by compression of the bounded diagonal and prime operator. For the residual e, a polynomial inverse candidate w_n=p_(n-1)(C)e has exact residual e_n=q_n(C)e, with q_n(x)=1-x p_(n-1)(x). Choose

    q_n(x)=T_n((u+ell-2x)/(u-ell))/T_n((u+ell)/(u-ell)).

The Chebyshev polynomials T_n are defined by T_0=1,T_1=x,T_(n+1)=2x T_n-T_(n-1). q_n(0)=1, so p_(n-1) is a polynomial. On [ell,u], |T_n|<=1. For z>1, T_n(z)=(s^n+s^-n)/2, s=z+sqrt(z^2-1), by the same recurrence. Therefore

    ||e_n|| <= 2 ((sqrt(u/ell)-1)/(sqrt(u/ell)+1))^n ||e||.

The identical bound holds in the C^-1 energy norm since the polynomial commutes with C^-1. Exact rational comparison u/ell < (53/10)^2 yields the simpler factor 2(43/63)^n. Upper bounds from direct integer powers:

[
  {
    "degree": 16,
    "decimal_ceiling": "1e-2"
  },
  {
    "degree": 32,
    "decimal_ceiling": "1e-5"
  },
  {
    "degree": 48,
    "decimal_ceiling": "1e-7"
  },
  {
    "degree": 64,
    "decimal_ceiling": "1e-10"
  },
  {
    "degree": 96,
    "decimal_ceiling": "1e-15"
  }
]

This is a bound for exact operator application only. A discretized or interval approximation cannot simply inherit it. If polynomial recurrence terms are approximated, use the recurrence's actual defects or certify the final e-Cw directly; the latter is already supported by the centered coefficient receiver. A separate higher-precision physical action of C on arbitrary piecewise functions would enable scalable exact-residual proposals without a dense full joint solve. Do not reuse the historical 5e-9 absolute error.

## 3. Scope

This reduces the coupled penalty to direction-specific signed data and a smaller auxiliary inverse error if a suitable w exists. It does not establish that the resulting finite correction or pole energy is favorable. No new runner activated, no actual polynomial candidate evaluated, no complete pole sign or global positivity conclusion.
