"""RETURN097-B: independently evaluate the first 400 mixed rows.
Use RETURN097-A's alternate shared scalars and interval coefficient
contraction, but compare with the distinct historical early-source bank.
The parent loader is used ONLY to read/rephase/normalize the reference.
No source or terminal acceptance on mere overlap.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
p=Path(__file__).resolve().parent.parent/'master-return-097a/regenerate_grouped.py'
source=p.read_text()
anchor='# Actual grouped tables; earlier 400 rows used different evaluation.'
assert source.count(anchor)==1
# Load just the alternate scalar and chart definitions, not the census.
ns={'__file__':str(p),'__name__':'alternate_grouped'}
exec(compile(source.split(anchor)[0],str(p),'exec'),ns)
w=ns['w'];S=ns['S']
ks=sorted(k for k,v in w['positive_even_weights'] if k>116)[:400]
assert len(ks)==400 and len(set(ks))==400
# First independently generate both alternate banks, before loading
# any historical mixed row code.
new=[[],[]]
for i,k in enumerate(ks):
 fac,ci,pi=ns['shared'](k)
 for par in (0,1):
  a=ns['calc'](par,k,fac[par],ci,pi)
  assert a.shape==(87,2) and np.all(a[:,0]<=a[:,1])
  new[par].append(a.tolist())
 if (i+1)%64==0:print('alternate early rows',i+1,k,flush=True)
(out/'early_regenerated.json').write_text(json.dumps({'ks':ks,'scale':S,'blocks':new},separators=(',',':'))+'\n')
# Historical reader is the same reference assembly used by RETURN083.
# Reset interval precision for its original context if it imports it.
sys.path.insert(0,str(root/'return067'))
import return067_contract as prev
assert list(prev.ks)==ks,('historical labels',len(prev.ks),len(ks))
records=[];fails=[]
for par in (0,1):
 old=np.asarray(prev.rawrows(par,S),dtype=object);a=np.asarray(new[par],dtype=object)
 assert old.shape==a.shape==(400,87,2) and np.all(old[:,:,0]<=old[:,:,1])
 contain=(old[:,:,0]<=a[:,:,0])&(a[:,:,1]<=old[:,:,1])
 over=(old[:,:,0]<=a[:,:,1])&(a[:,:,0]<=old[:,:,1])
 for i,j in np.argwhere(~contain)[:10]:
  fails.append({'parity':par,'k':ks[i],'retained_index':int(j),'old':old[i,j].tolist(),'new':a[i,j].tolist(),'overlap':bool(over[i,j])})
 rec={'parity':par,'rows':400,'entries':int(contain.size),'contained':int(contain.sum()),'overlap':int(over.sum())}
 records.append(rec);print('early mixed row comparison',rec,flush=True)
ok=all(r['contained']==r['entries'] for r in records)
(out/'receipt.json').write_text(json.dumps({'status':'PASS' if ok else 'NOT_CERTIFIED','scope':'first 400 mixed rows, alternate grouped scalars and frozen raw coefficient intervals; historical reader used for reference only','records':records,'mismatch_samples':fails,'alternate_evaluator_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2)+'\n')
assert ok,('containment not complete',fails[:3])
print('RETURN097-B FIRST 400 MIXED ROWS PASS',flush=True)
