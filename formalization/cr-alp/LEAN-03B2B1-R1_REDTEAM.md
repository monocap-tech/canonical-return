# LEAN-03B2B1-R1 pre-activation red-team

R0 native receipt: 0fccb454040dc47974d9f00b0dbc93178f6860e2, receipt blob 3e623d75903e505014d268dad3417698846f0ca6, root compiler log blob 6d5d2f5b941cc2cc66bb2cfca97df0d79f4bdcb3. Original Python and independent JS replay matched all 1600 terms, all endpoints and 3200 invalid endpoint rejection checks (replay blob fe0e0321999a12bc23198cf5d2dc388eec5b1876).

Root goal: integer cast of power/product versus power/product of real casts, in expTerm. Add Int.cast_mul and Int.cast_pow to the existing simp set. Both exact rules already passed in atanTerm in the same module and toolchain. No new theorem or assumption; expTerm statement unchanged. The remaining 5 downstream modules did not reach proof elaboration because RationalTrace.olean was absent; do not infer additional mathematical failures or successes.

Pre-submission checks executed by exact source comparison: unique proof site, all seven named audit directives unchanged, statement unchanged, forbidden proof constructs absent, both added cast rules present in passing atanTerm. All pass. Reversal of the patch recovers the exact R0 source. New source eaf14104c684d74663c0cc7c918ef03156aefab0 blob 3d448b5064c912afc2c0766fac9b5d086bda54d0. Generator and independent replay unchanged.

Runner must distinguish dependency-blocked modules from failed elaborations and require all 26 named audits, six compiled artifacts, exact original replay and durable remote custody before acceptance. No local native compilation available; static review is not kernel evidence. No change in mathematical scope.
