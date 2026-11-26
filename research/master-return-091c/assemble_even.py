"""RETURN091-C: independently assemble even TT/U and replay both frozen
RETURN086 terminal variants. Auxiliary saddle and raw H are shared inputs.
"""
import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve(); banks=Path(sys.argv[2]).resolve(); out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
O=10**12
def load(p):return json.loads(p.read_bytes())
def arr(x):return np.asarray(x,dtype=object)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scale(a,f):
 f=F(f);assert f>0
 l=a[:,:,0]*f.numerator//f.denominator
 h=-((-a[:,:,1]*f.numerator)//f.denominator)
 assert np.all(l<=h)
 return np.stack([l,h],axis=-1)
def compare(a,b,label,contain):
 assert a.shape==b.shape
 ok=(a[:,:,0]>=b[:,:,0])&(a[:,:,1]<=b[:,:,1]) if contain else (a[:,:,0]<=b[:,:,1])&(b[:,:,0]<=a[:,:,1])
 print(label,int(ok.sum()),ok.size,flush=True)
 assert np.all(ok),label
r87=root/'return087';r86=root/'return086'
pins={'finite_bank.json':'c6a1df49ff86806810e8763c10223b3b4b2e331b0243693802f1372eb14c443f','band_0.json':'624d811d17b570bdc58d39d7e47350bad0f91524443e8e2fad0cba91a7c6db58'}
for f,h in pins.items():assert sha(r87/f)==h
fb=load(r87/'finite_bank.json');bb=load(r87/'band_0.json')
assert fb['scale']==O and bb['record']['parity']==0
even=[v for v in fb['banks'] if v['parity']==0]
assert [(v['start'],v['end']) for v in even]==[(0,1680),(1680,2192),(2192,8111)]
p79=load(root/'return079/RETURN079_contract.json')
p81=load(root/'return081/contract.json');p80=load(root/'return080/band.json')
assert p79['scale']==p81['scale']==p80['scale']==O
assert p79['ks']==list(range(3520,4544,2))
assert p81['ks']==list(range(4544,16386,2))
assert p80['ks']==list(range(16386,32770,2))
high=load(banks/'089'/'bank_0.json');tail=load(out/'even_tail.json')
assert tail['K']==262144
source={}
for k,shape in [('tt',(87,87,2)),('U',(518,87,2))]:
 e,l,b=map(arr,(even[1]['cumulative'][k],even[2]['increment'][k],bb['bank'][k]))
 for x,p,label in [(e,p79['blocks'][0]['finite'][k],'2192 actual'),(l,p81['blocks'][0][k],'5919 actual'),(b,p80['blocks'][0][k],'8192 feature')]:
  assert x.shape==shape
  compare(x,arr(p),'even '+label+' '+k,True)
 h=scale(arr(high['matrices'][k]),F(O,high['scale']))
 t=scale(arr(tail['matrices'][k]),F(O,tail['scale']))
 assert h.shape==t.shape==shape
 source[k]=e+l+b+h+t
# Use the exact frozen RETURN086 verifier, including modular rank checks,
# exact doubled-center residual, raw H lower, and fixed congruence.
sys.path.insert(0,str(r86))
from check_terminal import check
required={'O','T','Q','XS','HS','CS','kl','ku','X','Y','V','C','Hraw','tt','U','target'}
captures=[]
for p in sorted(r86.glob('*.json')):
 d=load(p)
 if isinstance(d,dict) and required<=d.keys() and len(d['kl'])==518:
  assert F(d['target'])==F(49,50),(p.name,d['target'])
  captures.append((p,d))
assert len(captures)==2,[(p.name,d['target']) for p,d in captures]
records=[]
for p,d in captures:
 assert d['O']==O
 for k,m in source.items():
  compare(m,arr(d[k]),p.name+' original source overlap '+k,False)
  d[k]=m.tolist()
 q=out/('reconstructed_'+p.name);q.write_text(json.dumps(d,separators=(',',':')))
 rec=check(q);rec['parent_capture']=p.name;rec['parent_capture_sha256']=sha(p)
 records.append(rec)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'independent even source contractions; shared frozen auxiliary saddle and raw retained-energy source; RETURN086 exact checker','records':records,'pins':pins},indent=2)+'\n')
print('RETURN091-C BOTH EVEN TERMINAL COMPARISONS PASS',flush=True)
