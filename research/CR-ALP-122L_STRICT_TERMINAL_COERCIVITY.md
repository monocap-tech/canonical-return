# CR-ALP-122L — Strict terminal separation and conditional collar coercivity

Research-only. Independent exact integer shifted Sylvester tests on the four frozen 121 terminal lower-grid matrices at scale 10^18 (scale from 075_shifted_energy.py; 121 grid lower bounds K T K-K). Source commit 7c9dc0301df765ab31f9332e9e8c76015b59894c. No numerical optimization and no new physical integration.

[
  {
    "name": "odd",
    "dimension": 43,
    "hash": "ef6a2ce2b37654f5785c37d4d793030e843eeb075a1d0380bf6242cc283aa402",
    "shifts": [
      {
        "exp": 3,
        "pass": false
      },
      {
        "exp": 4,
        "pass": false
      },
      {
        "exp": 5,
        "pass": true
      }
    ],
    "certified_shift": "1e-5"
  },
  {
    "name": "second_batch_scalar",
    "dimension": 42,
    "hash": "0519c83b32982a3a2239f750b235ed3191428fe39c30ad6d21a54526d3bf672b",
    "shifts": [
      {
        "exp": 3,
        "pass": false
      },
      {
        "exp": 4,
        "pass": true
      }
    ],
    "certified_shift": "1e-4"
  },
  {
    "name": "second_batch_combined",
    "dimension": 42,
    "hash": "42c495d86072e538250e00069268a726555ad3fb17c803d5db7651904da32c3d",
    "shifts": [
      {
        "exp": 3,
        "pass": false
      },
      {
        "exp": 4,
        "pass": true
      }
    ],
    "certified_shift": "1e-4"
  },
  {
    "name": "second_batch_weighted",
    "dimension": 42,
    "hash": "07630a6c30fb74f1f180333065b6dcd19d96e2d7f0bca4016196877ef819e71b",
    "shifts": [
      {
        "exp": 3,
        "pass": false
      },
      {
        "exp": 4,
        "pass": false
      },
      {
        "exp": 5,
        "pass": true
      }
    ],
    "certified_shift": "1e-5"
  }
]

Algorithm: for each integer matrix N, subtract 10^(18-e) from each diagonal, then fraction-free symmetric elimination with exact divisibility checks and all strictly positive leading pivots. Scan e=3,...,12 until a passing shift. Thus N/10^18 >= delta I for the recorded delta. Failed trial shifts are inconclusive for the true unrounded matrix; do not infer a physical negative direction. Results here were executed in JS BigInt against the reopened frozen matrices.

## Transfer of strictness

Use the 122 generic simultaneous compression argument: B>=0 and U_l* B U_l<=A_l with A_l>0 imply B<=T^-1 if T=sum w_l U_l A_l^-1 U_l* is SPD and weights sum to one. For K>0, KTK-K>=delta I gives T>=K^-1+delta K^-2 >=(1+delta/kappa)K^-1 where kappa>=||K||>0. Inverting yields B<=(1+delta/kappa)^-1 K. Consequently K_C-B >= K-B >= delta/(kappa+delta) K>0 whenever K_C>=K. This strengthens the conditional conclusion to a positive quantitative separation in the original common finite coordinates. No assumption of invertible B.

For p=V z in the physical finite subspace and q in Q85 intersect D0, write M z=Xz+C^-1 r_C z, and c=inf spec C>0. The completed square gives W(p+q)>=m||z||^2+c||q+Mz||^2, where m=delta lambda_min(K)/(kappa+delta)>0. Since ||p+q||^2=||Vz||^2+||q||^2 <=(||V||^2+2||M||^2)||z||^2+2||q+Mz||^2, the full fixed-collar form has conditional coercivity W(f)>=min(m/(||V||^2+2||M||^2),c/2)||f||^2. The symbol V here denotes the actual physical coordinate map from finite Euclidean coordinates to E85, not just its coefficient matrix. This applies separately to each parity; take the lesser constant. It is a full infinite-dimensional coercivity implication, not a bound only on E85, and does not require Mz to lie in the log form domain.

The 122J geometry census and 122K source/trial/radius census have passed. Remaining assumptions: all original physical integrals genuinely enclosed under the stated midpoint arithmetic and higher-moment wrappers; source and cap comparators attach to the same signed form, and exact transport to the intended current original source. No numeric value for the full physical coercivity constant is claimed (requires rigorous kappa, lambda_min(K), ||V|| and ||M|| bounds). No fixed-collar promotion or arbitrary-support globalization.
