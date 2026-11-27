"""RETURN104-B exact-integer auxiliary Gram backend.
Transform only the finite dot-product backend of reviewed 092B/C.
All row generation, interval radii, tail bounds, signed metric and
terminal arithmetic stay identical. No binary64 dot products.
"""
import sys,json,hashlib,subprocess,zipfile
from pathlib import Path
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();work=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve()
work.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
base=Path(__file__).resolve().parents[1]
ref='87c9cd17b29da9d0a57cd426e0aa76a24eea0f65';prefix='archive/master-return-088-099'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
def get(p):return subprocess.check_output(['git','show',ref+':'+prefix+'/'+p])
assert subprocess.check_output(['git','rev-parse',ref+':'+prefix+'/manifest.json'],text=True).strip()=='37eb603b84a5a92b78d6623ab74fe880e191ef6c'
d=json.loads(get('manifest.json'))
labels={15038759268:'088',15142039146:'091c',15142054118:'092a',15148420383:'092b',15148031871:'092c'}
for a in d['archives']:
 if a['id'] not in labels:continue
 dest=work/labels[a['id']];dest.mkdir(exist_ok=True);z=dest/'original.zip'
 with z.open('wb') as f:
  for p in a['parts']:
   b=get(p['path']);assert len(b)==p['bytes'] and hashlib.sha256(b).hexdigest()==p['sha256'];f.write(b)
 assert z.stat().st_size==a['bytes'] and sha(z)==a['sha256']
 with zipfile.ZipFile(z) as f:
  assert f.testzip() is None;f.extractall(dest)
 z.unlink()
 print('authenticated',labels[a['id']],flush=True)
assert all((work/x).is_dir() for x in labels.values())
def replace(s,old,new):
 assert s.count(old)==1,old
 return s.replace(old,new)
def transform(path):
 s=path.read_text()
 s=replace(s,'import numpy as np','import numpy as np\nfrom flint import fmpz_mat')
 # Erase the original floating-point accumulation, keeping each
 # weighted midpoint and its directed interval radius unchanged.
 s=replace(s,' accum={(i,j):np.zeros((n,n),dtype=np.float64) for i in range(L) for j in range(i,L)}',' M_exact=fmpz_mat(n,n)')
 start=s.index('  ds=digits(mid,L)');end=s.index('  if start%4096',start)
 s=s[:start]+'''  block=fmpz_mat(mid.shape[0],mid.shape[1],[int(v) for v in mid.flat])
  M_exact += block.transpose()*block
'''+s[end:]
 start=s.index(' M=np.zeros((n,n),object)');end=s.index(' E=ev[:,None]',start)
 s=s[:start]+''' M=np.asarray([[int(M_exact[i,j]) for j in range(n)] for i in range(n)],dtype=object)
'''+s[end:]
 # Remove now-unused binary64 limb helper altogether.
 start=s.index('def digits(x,n):');end=s.index('\ndef ',start+1)
 # 092-C has an intervening bank path / proposal setup after
 # digits. Only delete through the literal return statement.
 end=s.index(' return ds\n',start)+len(' return ds\n')
 s=s[:start]+s[end:]
 # Keep the unchanged original interval and comparison algorithm.
 return s
def run(label,src,args):
 path=out/(label+'_exact.py');path.write_text(transform(src))
 subprocess.run([sys.executable,str(path),*map(str,args)],check=True)
 return {'source_sha256':sha(src),'transformed_sha256':sha(path)}
newb=out/'092b';newc=out/'092c';newb.mkdir();newc.mkdir()
r1=run('092b',base/'master-return-092b/common_saddle.py',[root,newb])
# Give augmented reconstruction the freshly exact common Gram, not the
# original common archive. Preserve other authenticated source banks.
import shutil
shutil.move(str(work/'092b'),str(work/'092b_original'))
shutil.copytree(newb,work/'092b')
r2=run('092c',base/'master-return-092c/augment_and_check.py',[root,newc,work])
records=[]
for p in (0,1):
 name=f'common_{p}.json';x=load(newb/name);y=load(work/'092b_original'/name)
 assert x['scale']==y['scale'] and x['nodes']==y['nodes'] and x['K']==y['K']
 n=518+p;records.append({'bank':name,'entries':n*n,'exact_equality':True})
x=load(newc/'odd_augmented.json');y=load(work/'092c'/'odd_augmented.json')
assert x['scale']==y['scale'] and x['K']==y['K']
records.append({'bank':'odd_augmented','entries':527**2,'exact_equality':True})
for p in (0,1):
 for suffix in ('','_alt'):
  name=f'reconstructed_matrices_{p}{suffix}.json'
  x=load(newc/name);y=load(work/'092c'/name)
  assert x==y,(p,suffix)
  print('exact backend terminal input equality',p,suffix or 'generator',flush=True)
assert load(newb/'receipt.json')['status']==load(newc/'receipt.json')['status']=='PASS'
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'exact FLINT integer finite auxiliary Gram products; identical original interval radii, tail formulas, final K and four terminal inputs','records':records,'source_transformations':[r1,r2],'exact_matrix_entries':sum(r['entries']for r in records)},indent=2)+'\n')
print('RETURN104-B EXACT INTEGER AUXILIARY GRAM AND TERMINAL EQUALITY PASS',flush=True)
