# CR-ALP-123AF — Frozen scalar amplification obstruction
W02 delegated continuation, same CR-ALP pole restoration jurisdiction. Predecessor 123AE execution 45951947443 mathematically completed, both SUFFICIENT_BOUND_INCONCLUSIVE; exact original archive recovered by custody-only run 45952736423, commit 630672f35e628aa03f3e7fd7a271cb601b1d4168, recovery.json blob 53887fac4f7e2a0288d4ca21ca976e865ccd3bdd. All eight remote chunk objects present; runner reconstructed ZIP from remote Git readback, verified archive and member hashes. Archive 130454707 bytes, SHA256 7b8fecd5b134b74ab83653d6d073a7cb4fcefe9d66552dc88fcc70f8caec75d7; result.json 280485282 bytes SHA256 9db87986b650df49caed887a7a1b47f864285979e7451e57e339c416e08225a0. No mathematical rerun. Independent connector read of recovery.json and eight chunk metadata; not an independent local reassembly of all 130 MB.

## Exact lower obstruction for the frozen scalar receiver
In 123AB notation let C=Q(Dcap+P)Q on the actual pole-null high space and Z=R-C Mhat with the existing E850 responses. By the cap and inherited original signed prime norm bound (122E prime.log blob d960e41188bbe6a53a26f80126dee8d7cc4cfa6d), 0<cI<=C<=MI, M=4.447938456+21/8=7.072938456. The cap is the 10^-9 downward rounded 2049 band of 068 (blob 1a7fe8629daa8b8e6cec185a79b6746d6b28975d), as used by 123AC/AE. The sign of P does not affect the norm upper bound. Thus C^-1>=I/M on Q. All physical and common-coordinate assumptions of 122R, 123AB and 123AC retained.

The actual receiver uses Dbar=L Delta L*>=D=K^-1+N^-1, with exact frozen LDL decomposition. Let A=Z L sqrt(Delta), and eta_bar=||A* C^-1 A||. This is the coupling quantity for the conservative frozen Dbar metric (equal in operator norm to Dbar^(1/2) Z* C^-1 Z Dbar^(1/2)); it is NOT a lower bound on amplification with a different metric or actual physical W inverse. For each standard column e_i, Parseval and C^-1>=I/M give
    eta_bar >= (Delta_i/M) sum_{j>85 matching parity}|(ZL)_{ji}|^2
            >= (Delta_i/(M G^2)) sum_{85<j<=2048 matching parity} max(0,|m_ji|-r_ji)^2.
The last bound uses only stored exact integer midpoint/radius pairs of the already Q-projected defect with grid G=10^60; no tail estimate, no diagonal inverse majorant, no unprojected representative, no floating eigenvalue. It is a lower bound on the actual capped-comparator coupling in the frozen metric, conditional on the prefix intervals' validity. Column selection is a diagnostic, not an optimization of the response bank.

For the existing sufficient receiver sigma>=E-U-(sqrt(b2)+sqrt(eta U))^2, using U>0 fixed, b2>=0, a necessary condition for this *certificate* to be positive is eta<E/U-1. Any valid eta upper bound for eta_bar exceeds the following exact certified lower floors:
| sector | zero-based column | eta_bar lower (decimal display) | certified rational floor | E/U-1 (display) | (1+lower)U/E (display) |
|---|---:|---:|---:|---:|---:|
| cosh | 0 | 206.4483737861353 | 206 | 1.3475511107722828 | 88.36799030027944 |
| sinh | 2 | 4.914725171504481 | 49/10 | 1.6471525681295625 | 2.2343726019855876 |
Both floors strictly exceed the exact necessary thresholds by rational comparison. Even perfect upper enclosure of eta_bar cannot make the frozen scalar receiver positive with the same U,E,b2. This is a method-local impossibility for the fixed conservative metric, response bank, and scalar U budget; not a sign or lower bound on the true pole loss, not an obstruction to different metrics, directional residual bounds, smaller U, different trials, or positivity itself.

## Reproducible exact check
Executed locally against the original 123AC artifact downloaded during W02-ALP-01. The exact 123AC result also exists at the current repository path below. No output generated that requires a large file; this executable source and results suffice to reproduce. Run with the raw result path as argument.
```python
import hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
b=Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(b).hexdigest()=='a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473'
d=json.loads(b)
assert d['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER'
assert set(d['sectors'])=={'1','2'}
M=F('4.447938456')+F(21,8)
for k,s in d['sectors'].items():
    n=44-int(k); G=int(s['prefix_grid'])
    assert G==10**60
    p=list(map(F,s['LDL_pivots']))
    assert len(p)==n and min(p)>0
    assert [j for j,_ in s['prefix_mid_radius']]==list(range(87 if k=='1' else 86,2049,2))
    assert all(len(row)==n and all(type(m)==int and type(r)==int and r>=0 for m,r in row) for j,row in s['prefix_mid_radius'])
    lo=[v*F(sum(max(0,abs(row[i][0])-row[i][1])**2 for j,row in s['prefix_mid_radius']),G*G)/M for i,v in enumerate(p)]
    i=max(range(n),key=lo.__getitem__)
    E,U=F(s['E']),F(s['U'])
    assert E>0 and U>0 and F(s['finite_charge_upper'])>=0
    floor=F(206) if k=='1' else F(49,10)
    assert i==(0 if k=='1' else 2)
    assert lo[i]>=floor>E/U-1
    print(k,i,float(lo[i]),str(floor),float(E/U-1),float((1+lo[i])*U/E),'EXACT OBSTRUCTION PASS')
```

## Reassessment and bounded successor
123AE diagnostics from recovered receipt: prefix spectral upper 1934.330247/111.298294, prefix squared box errors 4.2183e-56/4.1244e-56; tail leading spectral upper 1002.953936/62.964975, tail remainder Frobenius squared upper 3251.267549/200.217313; joint perturbation reserve 9431.716464/573.588561; best eta upper 9638.131033/599.039511. Upper-bound decomposition alone is not proof of intrinsic coupling. The new explicit prefix lower test, unlike the upper diagnostics, rules out polishing the frozen Dbar scalar certificate alone. STOP this method at the stated fixed inputs, not the pole restoration objective. No bank expansion.
An eligible next operation is assess direction-specific control of the actual vector Z* C^-1 e in the same metric, instead of sup over arbitrary residual directions. From 123AB exact identity a-M*r=b-Z*C^-1 e, a certified bound on that specific vector's Dbar norm may replace sqrt(eta U). Reuse the existing frozen trial, signed residual source and existing bank; first establish which signed pairings and independent error certificates are actually available. Do not assume a small pairing, or infer it from the upper-bound obstruction. If direction-specific evidence is insufficient, reevaluate the comparator metric and U budget rather than blindly enlarge the bank.

Four axes: mathematical conditional 122P/R unchanged, both pole signs OPEN; computational 123AE inconclusive and 123AF exact rational lower obstruction PASS; custodial original 123AE output reconstructed from remote chunks, 123AF reproducible text receipt; operational ordinary delegated continuation, no CI or new response solves for this check. No global claim.
