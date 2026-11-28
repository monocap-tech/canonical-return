# CR-ALP LEAN-03B1-R2 — Pre-submission red-team

Scope: repair one failed coercion elaboration, no mathematical claim promotion. Review performed before workflow activation.

## Immutable inputs and failure

- R1 custody: 5ca3d0644414627918a089c8bd7b4b2f1ba50f09; receipt blob fadab467f2a881dffabb12ad68d469fdf473ff18; compiler log blob f3ee09154eb05944b1e1cbe0b328855ad55d488c.
- Pinned Mathlib 5ed2965256430c3649e86755f9576b54eca72435, Mathlib/Analysis/SpecialFunctions/Log/Deriv.lean blob 56c2032d5fbaf0982487b6eff864a4f04bea9fa3, Real.sum_range_sub_log_div_le.
- R1 log shows remaining denominator 2 * ↑i + ↑1 versus 2 * ↑i + 1. The previous cast simp set omitted Nat.cast_one. Exponent 2*400+1 is separately normalized to 801 using decide.
- R2 source cb0f952be66e5a9753e63a6d7115d2281fd4df8e blob b9bf86656be650f89ffb39dc729f26491a20c220. TrigSeries unchanged (b8b092604b260ea0b0e0a36aaded778ecd7aa625).

## Adversarial checks (executed in JS BigInt and source comparison)

1. Unique proof site replaced; reversing replacement gives exactly R1 bytes: PASS.
2. All 15 elementary named audit directives unchanged: PASS.
3. No forbidden proof constructs in revised source: PASS.
4. Logarithm polynomial definition unchanged: PASS.
5. All declarations after logReduced byte-identical: PASS.
6. All 400 denominators are positive exact odd integers 1 through 799: PASS.
7. Last retained power 799, first omitted power 801: PASS.
8. For |u| <= 7/20, 1-u^2 >= 351/400 >= 1/2: PASS.
9. 4*(7/20)^801 < 10^-340: PASS.
10. Even 10 times that looser remainder is < 10^-340: PASS.

The exact arithmetic checks can be reproduced with Node:
```js
const assert=require('node:assert/strict');
for(let i=0;i<400;i++) assert(2n*BigInt(i)+1n===BigInt(2*i+1) && 2*i+1>0);
assert(2n*399n+1n===799n && 2n*400n+1n===801n);
assert(2n*(400n-49n)>=400n);
assert(4n*7n**801n*10n**340n<20n**801n);
assert(40n*7n**801n*10n**340n<20n**801n);
```

The review deliberately distinguishes mathematical guard sufficiency from elaboration success. Native Lean/Lake is unavailable in the local review environment; the above is a static source and independent arithmetic review, not a successful compilation or a substitute for native checking. The patched proof still requires actual kernel acceptance. R1's 19 passing named audits do not imply a closed elementary module. The new run must certify all 21 named declarations, emit both oleans, reject nonstandard axioms, and preserve the complete receipt remotely before 03B2. No change to numerical inputs, precision, or scope.
