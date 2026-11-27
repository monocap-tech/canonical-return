# CR-ALP-123K — Matrix-order inverse response improvement

Research-only. While 123I run 45922720120 executes, independently reopened 121 terminal result at immutable commit 7c9dc0301df765ab31f9332e9e8c76015b59894c (blob 2b55baf0e81ac9a296276f396f02633553e58b28). Use the actual terminal lower grid N, not merely the scalar m of 122P. All conclusions conditional on the common original comparator and simultaneous physical compression premises.

## Exact strengthened shifts

Fraction-free symmetric elimination with exact divisibility checks, all positive leading pivots, scale 10^18:
- Odd 43: N >= 10^-5 I; next tested 1.1e-5 is inconclusive.
- Even combined 42: N >= 2.1e-4 I; next tested 2.2e-4 is inconclusive.
Original matrix hashes respectively ef6a2ce2b37654f5785c37d4d793030e843eeb075a1d0380bf6242cc283aa402 and 42c495d86072e538250e00069268a726555ad3fb17c803d5db7651904da32c3d. The binary search and final neighbor shifts were checked directly in JS BigInt. The failed shifts concern the lower grid and are not negative physical witnesses.

## Preserve the matrix inverse geometry

From KTK-K >= N >0, T >= K^-1+K^-1 N K^-1 and B_C <= T^-1 <= G0=K(K+N)^-1 K. Therefore S_C=K_C-B_C >= S0=K-G0=(K^-1+N^-1)^-1. The last identity holds without commutativity: K-K(K+N)^-1 K=K(K+N)^-1 N, whose inverse is N^-1(K+N)K^-1=N^-1+K^-1.

The matrix-order receiving bound is stronger than separately charging kappa/m. Put A=K^-1/2 N K^-1/2. Then congruence by K^-1/2 sends G0 to (I+A)^-1 and S0 to A(I+A)^-1. Hence if N>=delta I and K<=kappa I, A>=delta/kappa I implies G0<=(kappa/delta) S0. This is an exact Loewner implication, not an eigenvalue approximation or an identity between noncommuting matrices.

For the 123F auxiliary residual e, U>=<e,C^-1 e>, and b_y=a-X*r-r_C*y as defined there, the transformed finite response norm satisfies

    ||S0^-1/2 r_C* C^-1 e||^2 <= (kappa/delta) U.

Proof: the map S0^-1/2 r_C* C^-1/2 has squared norm bounded by that of S0^-1/2 B_C S0^-1/2 <= kappa/delta I. Also S0^-1=K^-1+N^-1<= (1/lambda+1/delta) I. Consequently the actual pole Schur scalar has lower bound

    joint_energy - U - (sqrt((1/lambda+1/delta))*||b_y|| + sqrt((kappa/delta)*U))^2.

The previous 123I loss was U+(||b_y||+sqrt(kappa U))^2/m. Both squared terms now have smaller sufficient scalar coefficients. This does not require recovering T or inverting the physical C.

Independent exact rational ceilings from 122P:
[
  {
    "k": "1",
    "delta": "1/100000",
    "finite_old": "1986582",
    "finite_new": "100010",
    "error_old": "3946485",
    "error_new": "198658"
  },
  {
    "k": "2",
    "delta": "21/100000",
    "finite_old": "212009",
    "finite_new": "4772",
    "error_old": "449456",
    "error_new": "10096"
  }
]

The stored 123I output reports U, ||b_y|| upper and joint energy lower. Therefore this stronger bound can be checked from those receipts without any new special-function evaluation, even if 123I's original sufficient test is inconclusive. An even stronger later receiver can compute b_y* (K^-1+N^-1) b_y directly with certified coordinates and determine the generalized ratio N>=t K instead of using scalar extrema.

## Boundary

No new runner or unconditional proof. These are exact arithmetic shift checks and a classical conditional Loewner transfer. A positive complete pole lower bound requires actual 123I error enclosures; no sign is claimed. Global support and RH remain open.
