"""RETURN092-A: independent odd terminal cross-check and exact source inventory.
No new trial solve. Do not reuse legacy odd terminal arithmetic or silently
promote the shared auxiliary saddle or retained-energy inputs.
"""
import sys,json,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve(); bank=Path(sys.argv[2]).resolve(); out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
def load(p):return json.loads(p.read_bytes())
def arr(x):return np.asarray(x,dtype=object)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# Record exact archived source, not snippets guessed from filenames.
# Include all Python files from the relevant auxiliary construction stages.
inventory={}
for dirname in ('return057','return062','return063','return077','return081','return083','return086','return087'):
 directory=root/dirname
 assert directory.is_dir(),dirname
 for p in sorted(directory.glob('*.py')):
  rel=str(p.relative_to(root));inventory[rel]={'sha256':sha(p),'bytes':p.stat().st_size}
  dest=out/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(p.read_bytes())
  if dirname in ('return063','return083') and ('inverse' in p.name or 'constraint' in p.name or p.name=='setup.py'):
   print('BEGIN AUTHENTICATED SOURCE',rel,sha(p),flush=True)
   print(p.read_text(),flush=True)
   print('END AUTHENTICATED SOURCE',rel,flush=True)
# Inspect actual serialized input field types/shapes without emitting whole matrices.
schemas={}
def describe(x,depth=0):
 if isinstance(x,dict):return {k:describe(v,depth+1) for k,v in x.items()} if depth<3 else {'keys':list(x)}
 if isinstance(x,list):return {'length':len(x),'first':describe(x[0],depth+1) if x and depth<5 else None}
 if isinstance(x,int):return x if abs(x)<=10**55 else 'large integer'
 if isinstance(x,str):return x if len(x)<100 else 'long string'
 return type(x).__name__
for dirname in ('return057','return062','return063','return077'):
 for p in sorted((root/dirname).glob('*.json')):
  d=load(p);rel=str(p.relative_to(root))
  schemas[rel]={'sha256':sha(p),'schema':describe(d)}
(out/'source_inventory.json').write_text(json.dumps({'sources':inventory,'schemas':schemas},indent=2)+'\n')
print('source inventory',len(inventory),'Python files',len(schemas),'JSON records',flush=True)
# Check actual archived odd independent assembly, not original aggregate.
a=load(bank/'assembled.json');receipt=load(bank/'receipt.json')
assert receipt['status']=='PASS' and F(receipt['target'])==F(19,20)
assert a['source_scale']==10**12
K=arr(a['K']);tt=arr(a['tt']);U=arr(a['U'])
assert K.shape==(527,527,2) and tt.shape==(87,87,2) and U.shape==(527,87,2)
for x in (K,tt,U):assert np.all(x[:,:,0]<=x[:,:,1])
r86=root/'return086';sys.path.insert(0,str(r86))
from check_terminal import check
records=[]
for name in ('matrices_1.json','matrices_1_alt.json'):
 p=r86/name;d=load(p)
 assert F(d['target'])==F(19,20) and d['O']==a['source_scale']
 for key,m in [('tt',tt),('U',U)]:
  old=arr(d[key]);assert old.shape==m.shape
  ok=(old[:,:,0]<=m[:,:,1])&(m[:,:,0]<=old[:,:,1])
  print(name,'source overlap',key,int(ok.sum()),ok.size,flush=True);assert np.all(ok)
 oldl,oldu=arr(d['kl']),arr(d['ku'])
 assert oldl.shape==oldu.shape==(527,527)
 # Compare before changing scale, exact integer products, no intersection.
 ok=(oldl*a['K_scale']<=K[:,:,1]*d['T'])&(K[:,:,0]*d['T']<=oldu*a['K_scale'])
 print(name,'auxiliary saddle overlap',int(ok.sum()),ok.size,flush=True);assert np.all(ok)
 d['T']=a['K_scale'];d['kl']=K[:,:,0].tolist();d['ku']=K[:,:,1].tolist()
 d['tt']=tt.tolist();d['U']=U.tolist()
 q=out/('reconstructed_'+name);q.write_text(json.dumps(d,separators=(',',':')))
 rec=check(q);rec['parent_capture_sha256']=sha(p);records.append(rec)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'independent odd source assembly and alternate terminal arithmetic; auxiliary baseline and raw retained-energy source shared','records':records,'input_sha256':sha(bank/'assembled.json'),'source_inventory_sha256':sha(out/'source_inventory.json')},indent=2)+'\n')
print('RETURN092-A BOTH ODD TERMINAL VARIANTS PASS',flush=True)
