"""RETURN105-A staged literal high-band anchor and scalar attachment.
Consume authenticated original files and the permanent RETURN101-A bank.
No new quadrature, terminal solve, or mathematical promotion.
"""
import sys,json,hashlib,zipfile,ast
from pathlib import Path
from fractions import Fraction as F
from flint import arb,acb,fmpq,ctx
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();archive=Path(sys.argv[2]).resolve()
bank=Path(sys.argv[3]).resolve();out=Path(sys.argv[4]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
def A(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
with zipfile.ZipFile(archive) as z:
 names=[n for n in z.namelist() if Path(n).name=='archive_manifest.json'];assert len(names)==1
 b=z.read(names[0])
assert hashlib.sha256(b).hexdigest()=='38b4781ddac2b4f0437f75f8eb802e78a54c1fbc366f4690880aa73b0aa9d112'
manifest=json.loads(b);assert len(manifest)==383
source=root/'return073/return073_band.py'
assert sha(source)==manifest[str(source.relative_to(root))]
# Authenticate all of the original reproduction files before importing
# reference initialization. RETURN086 direct_sums has its separate pin
# in the independently accepted RETURN101-A direct scalar receipt.
for name,digest in manifest.items():assert sha(root/name)==digest,name
c=load(bank/'constants.json');d=load(bank/'direct.json')
assert c['status']==d['status']=='PASS'
anchor=next(v for v in c['anchors'] if v['k']==32768)
old=tuple(F(x,10**55) for x in anchor['historical_units_1e55'])
ctx.prec=384
pi=arb.pi();R=2*pi*32768/5
g0=acb(A(F(1,4))).digamma().real-pi.log()
D0=acb(A(F(1,4)),R/2).digamma().real-pi.log()+g0+4
assert A(old[0])<D0<A(old[1])
# Reference import is solely to identify the literal scalar endpoints
# used by both independent high-band generating scripts.
sys.argv=[str(source)];sys.path.insert(0,str(source.parent))
import return073_band as ref
assert ref.K==32768 and ref.U==262144 and ref.Q==10**55
assert tuple(ref.D0)==old,('actual high-band anchor differs',str(ref.D0),str(old))
assert A(ref.pi[0])<pi<A(ref.pi[1])
assert ref.pi==(F(314159265358979323846,10**20),F(314159265358979323847,10**20))
base=Path(__file__).resolve().parents[1]
records=[]
for p in (base/'master-return-088/audit_high_band.py',base/'master-return-090/audit_restricted.py'):
 t=ast.parse(p.read_text())
 def val(name):
  a=[n.value for n in t.body if isinstance(n,ast.Assign) and any(isinstance(x,ast.Name) and x.id==name for x in n.targets)]
  assert len(a)==1,(p,name)
  return a[0]
 assert ast.dump(val('D0'))==ast.dump(ast.parse('ref.D0',mode='eval').body)
 assert eval(compile(ast.Expression(val('pi')),str(p),'eval'),{'F':F})==ref.pi
 assert ast.literal_eval(val('col'))==[[0,1,3,4],[1,2,5,6],[3,5,7,8],[4,6,8,9]]
 assert ast.literal_eval(val('K'))==32768
 records.append({'source':str(p.relative_to(base)),'sha256':sha(p),'D0_exact_equality':True,'pi_exact_equality':True,'scalar_product_index_identity':True})
def subset(x,s,y,t):
 assert len(x)==len(y)
 if isinstance(x[0],int):
  assert len(x)==len(y)==2 and y[0]*s<=x[0]*t<=x[1]*t<=y[1]*s
  return 1
 return sum(subset(a,s,b,t) for a,b in zip(x,y))
p=root/'return086/direct_sums.json'
assert sha(p)==d['historical_sha256']
oldbanks=load(p)['bands']
count=0
for v in d['bands']:
 b=next(x for x in oldbanks if x['name']==v['name'])
 assert v['p']==b['p']==list(range(2,47)) and v['count']==b['count']
 count+=subset(v['sums'],v['scale'],b['sums'],b['scale'])
assert count==900
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'literal D0 and pi identity in both high-band contractions; independent Arb anchor and direct scalar containment; no new solve','original_manifest_entries':383,'reference_source_sha256':sha(source),'D0':[str(x) for x in ref.D0],'pi':[str(x) for x in ref.pi],'scalar_containment':count,'records':records},indent=2)+'\n')
print('RETURN105-A LITERAL HIGH BAND SCALAR ATTACHMENT PASS',count,flush=True)
