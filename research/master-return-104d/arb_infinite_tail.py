"""RETURN104-D staged Arb infinite leading-tail scalar substitution.
No new analytic inequalities, contractions, or trial solutions. Re-evaluate
every literal phase chord and Gamma anchor in the accepted infinite-tail
source, preserving the exact rational contraction algorithm. Require exact
equality of all terminal-consumed interval banks to the custody originals.
"""
import sys,hashlib,json,subprocess,zipfile
from pathlib import Path
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
base=Path(__file__).resolve().parents[1]
ref='87c9cd17b29da9d0a57cd426e0aa76a24eea0f65';prefix='archive/master-return-088-099'
def get(p):return subprocess.check_output(['git','show',ref+':'+prefix+'/'+p])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
assert subprocess.check_output(['git','rev-parse',ref+':'+prefix+'/manifest.json'],text=True).strip()=='37eb603b84a5a92b78d6623ab74fe880e191ef6c'
manifest=json.loads(get('manifest.json'))
plans=[
 ('even',base/'master-return-091c/even_tail.py','08f85186c54b034fe89df628024ec09e5e3846b8',15142039146,'even_tail.json'),
 ('odd',base/'master-return-090b/audit_asymptotic.py','7f3c466e3905ffdeac2c9d81207183977476b7c5',15092621662,'asymptotic_bank.json')]
def replace(s,old,new):
 assert s.count(old)==1,old
 return s.replace(old,new)
records=[]
for label,path,blob,artifact,bank in plans:
 assert subprocess.check_output(['git','hash-object',str(path)],text=True).strip()==blob
 s=path.read_text()
 s=replace(s,'from mpmath import iv','from flint import arb,acb,fmpq,ctx')
 start=s.index('def ivf(x):');end=s.index('pi=pair(iv.pi)',start)
 s=s[:start]+'''ctx.prec=384
def arbq(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def interval(a):
 assert a[0]<=a[1]
 return arbq(a[0]).union(arbq(a[1]))
def pair(v):
 l=(v.lower()*S).floor().unique_fmpz()
 u=(v.upper()*S).ceil().unique_fmpz()
 assert l is not None and u is not None
 return F(int(l),S),F(int(u),S)
PI=arb.pi()
'''+s[end:]
 s=replace(s,'pi=pair(iv.pi)','pi=pair(PI)')
 s=replace(s,'iv.ln(p)/iv.sqrt(n)','arb(p).log()/arb(n).sqrt()')
 s=replace(s,'4*iv.pi*iv.ln(ivf(z))/5','4*PI*arbq(z).log()/5')
 s=replace(s,'2*abs(iv.sin(a/2))','2*abs((a/2).sin())')
 if label=='odd':
  start=s.index('sys.path[:0]=');end=s.index('graw=[]',start)
  s=s[:start]+'''R=2*PI*K/5
g0=acb(arbq(F(1,4))).digamma().real-PI.log()
gR=acb(arbq(F(1,4)),R/2).digamma().real-PI.log()
'''+s[end:]
  s=replace(s,'gR+src.G0+4','gR+g0+4')
  s=replace(s,'/(2*iv.pi)','/(2*PI)')
 assert 'iv.' not in s and 'mpmath' not in s and 'fast.sums' not in s
 target=out/label;target.mkdir()
 code=target/'arb_tail.py';code.write_text(s)
 subprocess.run([sys.executable,str(code),str(root),str(target)],check=True)
 new=load(target/bank)
 a=next(x for x in manifest['archives'] if x['id']==artifact)
 archive=target/'reference.zip'
 with archive.open('wb') as f:
  for part in a['parts']:
   b=get(part['path']);assert len(b)==part['bytes'] and hashlib.sha256(b).hexdigest()==part['sha256'];f.write(b)
 assert archive.stat().st_size==a['bytes'] and sha(archive)==a['sha256']
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  names=[n for n in z.namelist() if Path(n).name==bank];assert len(names)==1,names
  old=json.loads(z.read(names[0]))
 archive.unlink()
 # Comparison against the previously accepted infinite-tail arithmetic
 # enclosures, not against the much broader historical parent.
 fields=('scale','K','matrices') if label=='even' else ('scale','K','lambda','old','signed_reweighting','new_restoration_cross')
 for key in fields:
  assert new[key]==old[key],(label,key,'requires re-binding, not an automatic mathematical rejection')
 records.append({'parity':label,'reference_archive':artifact,'source_blob':blob,'transformed_sha256':sha(code),'bank_sha256':sha(target/bank),'equal_fields':list(fields)})
 print('Arb infinite tail exact interval equality',label,flush=True)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'Arb substitution for all literal infinite leading tail phase weights, chords and Gamma anchor; unchanged rational contractions and exact terminal-consumed bank equality','records':records},indent=2)+'\n')
print('RETURN104-D ARB INFINITE LEADING TAIL EQUALITY PASS',flush=True)
