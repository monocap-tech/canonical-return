"""RETURN093-A: reopen original retained-energy source and execute its
archived independent replay; bind terminal charts to RETURN059. This is
not yet a new scalar-level reconstruction of Hraw.
"""
import sys,json,hashlib,subprocess
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();bank=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
rawfile=root/'return059/RETURN059_raw_receipt.json'
assert sha(rawfile)=='2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8'
raw=load(rawfile)
# Preserve source text and exact digests, and disclose all relevant code
# before invoking the archived independent mathematical replay.
sources={}
for p in sorted(root.glob('return*/*.py')):
 s=p.read_text()
 if p.parent.name in ('return048','return059','return086') or ('Hraw' in s and p.parent.name in ('return057','return058','return064','return081','return083')):
  rel=str(p.relative_to(root));sources[rel]=sha(p)
  dst=out/'sources'/rel;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(p.read_bytes())
  if p.parent.name in ('return048','return059','return086') and ('raw' in p.name or p.name=='capture.py' or 'retained' in p.name):
   print('BEGIN RETAINED SOURCE',rel,sha(p),flush=True);print(s,flush=True);print('END RETAINED SOURCE',rel,flush=True)
print('retained source inventory',len(sources),flush=True)
candidates=[p for p in (root/'return059').glob('*.py') if 'PASS independent raw parity' in p.read_text() and 'PASS complete 174-dimensional original bank' in p.read_text()]
assert len(candidates)==1,[p.name for p in candidates]
replay=candidates[0];print('executing archived independent retained replay',replay.name,sha(replay),flush=True)
# The archive is extracted in isolated runner storage. Its original ZIP
# remains unchanged. Require the exact original source pins after replay.
run=subprocess.run([sys.executable,str(replay)],cwd=replay.parent,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
(out/'retained_replay.log').write_text(run.stdout)
print(run.stdout,flush=True);assert run.returncode==0
assert 'PASS complete 174-dimensional original bank' in run.stdout
assert sha(rawfile)=='2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8'
assert all(sha(root/p)==h for p,h in sources.items())
# Bind the exact rational chart, and confirm that RETURN092-C kept the
# literal RETURN086 raw energy array and scale in all four terminal variants.
# Hraw scalar derivation is separately audited in the next gate.
records=[]
for par in (0,1):
 for suffix in ('','_alt'):
  name=f'matrices_{par}{suffix}.json'
  a=load(root/'return086'/name);p=bank/('reconstructed_'+name);d=load(p)
  C=np.asarray(d['C'],object);c=np.asarray(raw['blocks'][par]['C'],object)
  assert C.shape==c.shape==(87,87)
  assert np.array_equal(C*raw['T'],c*d['CS'])
  assert d['HS']==a['HS'] and d['Hraw']==a['Hraw']
  assert np.array_equal(np.asarray(a['C'],object)*d['CS'],C*a['CS'])
  print('retained chart and frozen H custody',par,suffix,7569,7569,flush=True)
  records.append({'parity':par,'variant':suffix or 'generator','chart_entries':7569,'H_entries':7569,'H_scale':d['HS'],'terminal_sha256':sha(p),'capture_sha256':sha(root/'return086'/name)})
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'archived independent raw positivity replay and terminal chart/energy custody; no new scalar regeneration','source_pins':sources,'records':records,'raw_receipt_sha256':sha(rawfile)},indent=2)+'\n')
print('RETURN093-A RETAINED REPLAY AND TERMINAL CHART CUSTODY PASS',flush=True)
