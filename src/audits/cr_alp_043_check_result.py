#!/usr/bin/env python3
"""Small independent exact acceptance-ledger check (not coefficient reconstruction)."""
from fractions import Fraction as F
from pathlib import Path
import json
p=Path(__file__).resolve().parent
r=json.loads((p/'CR-ALP-043_exact_band_result.json').read_text())
t=json.loads((p/'CR-ALP-042_exact_result.json').read_text())
mu={x['N']:F(x['certified_margin']) for x in t['cutoffs']}
rows={(x['N'],x['parity']):x for x in r['results']}
assert set(rows)=={(N,k) for N in (172,242,300) for k in (1,2)}
for (N,k),x in rows.items():
 b,z,v=map(F,[x['band_upper'],x['tail_upper'],x['rho_upper']]);assert b>=0 and z>=0 and b+z==v
 assert F(x['mu'])==mu[N] and x['pass_positive']==(v<mu[N])
 assert x['band_rows']==7414 and x['band_end']==15000
 if N==172:assert v>=mu[N] # failure of this sufficient budget only
 else:
  eta={242:F(3,10),300:F(1,4)}[N];assert v<eta<mu[N]
  assert mu[N]-eta=={242:F(3,10),300:F(1,2)}[N]
print('PASS rational ledger: both parities at N=242,300; tail margins 3/10 and 1/2; N=172 inconclusive')
