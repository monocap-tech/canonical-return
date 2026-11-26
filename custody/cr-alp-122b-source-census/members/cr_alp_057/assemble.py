#!/usr/bin/env python3
"""Transport all saved coefficient radii to the complete projected physical residual.
The exact physical projector and prime shifts are retained; no moment comparison is evaluated.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,hashlib
ROOT=Path(__file__).resolve().parent
raw=(ROOT/'result.json').read_bytes();graw=(ROOT/'grouped_proxy.json').read_bytes();traw=(ROOT/'trial.json').read_bytes();a=json.loads(raw);g=json.loads(graw);x=json.loads(traw);S=int(a['scale']);assert S==int(g['scale'])==int(x['scale']);assert hashlib.sha256(raw).hexdigest()==g['source_result_sha256']
def radius(rows):
 z=0
 for row in rows:
  for l,h in row:
   l,h=int(l),int(h);m=(l+h)//2;z+=max(m-l,h-m)**2
 return F(isqrt(z)+1,S)
out={}
for k,d in a['sectors'].items():
 t=x['sectors'][k];dr=radius(d['source_sine_coefficients']);fr=radius([d['pole_coefficient']]);Dr=radius(t['D']);Dbr=radius(t['Db']);
 # Physical sine functions are orthonormal; ||b||<2 for both parities.
 # P norm <=21/8; K_{1/2} norm <=4 by the integral kernel Schur bound.
 src=F(21,8)*(dr+2*fr)+8*fr
 tr=F(1101,200)*(Dr+2*Dbr)
 # Use independent round-number Gamma bound: 4.480e-9 / 4.678e-9.
 gamma=F(4480 if k=='1' else 4678,10**12); total=gamma+src+tr;assert total<F(5,10**9)
 out[k]={'source_and_resonant_rounding_upper':str(src),'trial_rounding_upper':str(tr),'gamma_boundary_proxy_upper':str(gamma),'complete_residual_proxy_error_upper':str(total)}
 print(k,'source/resonance error',float(src),'trial rounding',float(tr),'complete projected residual error',float(total),flush=True)
(ROOT/'assembled_summary.json').write_text(json.dumps({'source_sha256':hashlib.sha256(raw).hexdigest(),'group_sha256':hashlib.sha256(graw).hexdigest(),'trial_sha256':hashlib.sha256(traw).hexdigest(),'sectors':out},indent=2)+'\n')
