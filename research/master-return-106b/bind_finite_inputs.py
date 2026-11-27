"""RETURN106-B: bind the actual RETURN087 consumed finite input rows to
the independent RETURN104-A Arb intervals, before reciprocal weighting.
No new solve; retain the verified exact integer contraction backend.
"""
import ast,hashlib,json,os,subprocess,sys,zipfile
from pathlib import Path
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
ref='4399c20b2b708260463dffa25f02797183d1983f';base='archive/master-return-103-104'
def get(p):return subprocess.check_output(['git','show',ref+':'+base+'/'+p])
assert subprocess.check_output(['git','rev-parse',ref+':'+base+'/manifest.json'],text=True).strip()=='47e3d9d70fdda55f4b0b2c76b1384101e394915d'
manifest=json.loads(get('manifest.json'))
bank=out/'arb';bank.mkdir()
a=next(x for x in manifest['artifacts'] if x['id']==15162478276)
z=out/'arb.zip'
with z.open('wb') as f:
 for p in a['parts']:
  b=get(p['path']);assert len(b)==p['bytes'] and hashlib.sha256(b).hexdigest()==p['sha256'];f.write(b)
assert z.stat().st_size==a['bytes'] and sha(z)==a['sha256']=='b1afca72214526a3b9581d0cbd546b5b783d672cdeb4d993be7f168b02824a84'
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None;f.extractall(bank)
z.unlink()
rec=load(bank/'receipt.json')
assert rec['status']=='PASS' and rec['entries']==rec['contained']==2001138
src=root/'return087/finite_audit.py'
assert sha(src)=='cf0f1170ecaf0d5faf8b4f1560517d8f678aeba9cb000f9baf68a42fd1b0c89c'
oldbank=src.parent/'finite_bank.json'
assert sha(oldbank)=='c6a1df49ff86806810e8763c10223b3b4b2e331b0243693802f1372eb14c443f'
original=load(oldbank)
s=src.read_text();print('BEGIN AUTHENTICATED FINITE SOURCE',flush=True);print(s,flush=True);print('END AUTHENTICATED FINITE SOURCE',flush=True)
tree=ast.parse(s)
func={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
assert 'dot' in func and 'weighted' in func
n=func['dot'];assert [v.arg for v in n.args.args]==['a','b']
assert [v.arg for v in func['weighted'].args.args][:2] in (['x','ids'],['a','ids'],['t','ids'],['x','ks'],['a','ks'],['t','ks']) or len(func['weighted'].args.args)==2
# Insert observational wrapper after the original weighted definition.
# The driver calls weighted with the actual unweighted input rows and
# ordered coordinate labels. Compare before any reciprocal multiplication.
wrapper='''
_return106_weighted_original=weighted
_return106_arb_path=Path(__import__('os').environ['RETURN106_ARB_PATH'])
_return106_receipt_path=Path(__import__('os').environ['RETURN106_RECEIPT_PATH'])
_return106_banks=[{},{}]
_return106_early=set()
for _p in sorted(_return106_arb_path.glob('*.json')):
 if _p.name in ('receipt.json','features.json'):continue
 _d=json.loads(_p.read_bytes())
 assert _d['scale']==10**12 and len(_d['blocks'])==2
 if _p.name=='early.json':
  assert len(_d['ks'])==400
  _return106_early=set(_d['ks'])
 for _par in (0,1):
  assert len(_d['ks'])==len(_d['blocks'][_par])
  for _k,_row in zip(_d['ks'],_d['blocks'][_par]):
   assert _k not in _return106_banks[_par] and len(_row)==87
   _return106_banks[_par][_k]=_row
assert len(_return106_early)==400 and all(len(b)==8111 for b in _return106_banks)
_return106_seen=[[],[]]
_return106_mismatch=[]
_return106_counts=[0,0]
def weighted(x,ids):
 # Determine the parity from the actual driver's loop, not a
 # favorable numerical comparison.
 par=sys._getframe(1).f_locals['par']
 assert par in (0,1) and S==10**16
 assert x.shape==(len(ids),87,2)
 for i,k in enumerate(ids):
  k=int(k);assert k in _return106_banks[par]
  _return106_seen[par].append(k)
  row=_return106_banks[par][k]
  for j,(l,u) in enumerate(row):
   a,b=map(int,x[i,j])
   ok=a*10**12<=l*S<=u*S<=b*10**12
   _return106_counts[par]+=int(ok)
   if not ok and len(_return106_mismatch)<30:
    _return106_mismatch.append([par,k,j,[a,b],[l,u]])
 return _return106_weighted_original(x,ids)
'''
exact='''def dot(a,b):
 from flint import fmpz_mat
 assert a.shape[0]==b.shape[0] and len(a)<=64
 A=fmpz_mat(a.shape[0],a.shape[1],[int(v) for v in a.flat])
 B=fmpz_mat(b.shape[0],b.shape[1],[int(v) for v in b.flat])
 C=A.transpose()*B
 return np.asarray([[int(C[i,j]) for j in range(b.shape[1])] for i in range(a.shape[1])],dtype=object)
'''
lines=s.splitlines(keepends=True)
edits=[(n.lineno-1,n.end_lineno,exact),(func['weighted'].end_lineno,func['weighted'].end_lineno,wrapper)]
for i,j,t in sorted(edits,reverse=True):lines[i:j]=[t+'\n']
end='''
assert all(len(v)==8111 and v==sorted(_return106_banks[p]) for p,v in enumerate(_return106_seen))
assert all(set(v[:400])==_return106_early for v in _return106_seen)
ok=all(v==8111*87 for v in _return106_counts)
_return106_receipt_path.write_text(json.dumps({'status':'PASS' if ok else 'NOT_CERTIFIED','scope':'literal consumed RETURN087 unweighted mixed rows against independent Arb bank at exact rational scales','rows_per_parity':8111,'early_rows_per_parity':400,'entries':8111*174,'contained':sum(_return106_counts),'by_parity':_return106_counts,'mismatch_samples':_return106_mismatch},indent=2)+'\\n')
assert ok,_return106_mismatch[:3]
print('RETURN106-B LITERAL FINITE SOURCE INPUT CONTAINMENT PASS',sum(_return106_counts),flush=True)
'''
new=''.join(lines)+'\n'+end
ast.parse(new)
instrumented=src.parent/'finite_audit_attached.py';instrumented.write_text(new)
(out/'instrumented_source.py').write_text(new)
env=dict(os.environ,RETURN106_ARB_PATH=str(bank),RETURN106_RECEIPT_PATH=str(out/'receipt.json'))
subprocess.run([sys.executable,str(instrumented),str(root)],cwd=src.parent,env=env,check=True)
assert load(oldbank)==original
r=load(out/'receipt.json');assert r['status']=='PASS'
r['original_source_sha256']=sha(src);r['instrumented_source_sha256']=sha(instrumented);r['original_finite_contraction_bank_exact_equality']=True;r['early_entries']=69600
(out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print('RETURN106-B EXACT FINITE CONTRACTION OUTPUT EQUALITY PASS',flush=True)
