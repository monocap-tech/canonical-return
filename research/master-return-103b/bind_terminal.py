"""RETURN103-B: literal terminal lineage from permanent custody.
No Actions artifact downloads, no new solve, no proof promotion.
"""
import hashlib,json,subprocess,sys,zipfile
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
work=Path(sys.argv[1]);out=Path(sys.argv[2]);work.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
def show(ref,path):return subprocess.check_output(['git','show',ref+':'+path])
refs=[
 ('87c9cd17b29da9d0a57cd426e0aa76a24eea0f65','archive/master-return-088-099','37eb603b84a5a92b78d6623ab74fe880e191ef6c','archives'),
 ('bce3d5491d0cabe684e05237ca3e0533f22a5914','archive/master-return-101-102','337ba07e5d56152627d1ccd3c19ca397f12740df','artifacts')]
labels={15142039146:'091c',15142054118:'092a',15148420383:'092b',15148031871:'092c',15149302173:'093b',15158236154:'099b',15158476209:'101a',15158421363:'101b',15158272662:'101c',15158975398:'102a',15158606500:'102b'}
banks=work/'banks';banks.mkdir(exist_ok=True)
authenticated=[]
for ref,base,blob,key in refs:
 assert subprocess.check_output(['git','rev-parse',ref+':'+base+'/manifest.json'],text=True).strip()==blob
 d=json.loads(show(ref,base+'/manifest.json'))
 assert len(d[key])==(20 if key=='archives' else 5)
 for a in d[key]:
  if a['id'] not in labels:continue
  name=labels[a['id']];z=work/(name+'.zip')
  with z.open('wb') as f:
   if key=='archives':
    for p in a['parts']:
     b=show(ref,base+'/'+p['path']);assert len(b)==p['bytes'] and hashlib.sha256(b).hexdigest()==p['sha256']
     f.write(b)
   else:f.write(show(ref,base+'/'+a['path']))
  assert z.stat().st_size==a['bytes'] and sha(z)==a['sha256']
  with zipfile.ZipFile(z) as f:
   assert f.testzip() is None
   f.extractall(banks/name)
  authenticated.append({'bank':name,'sha256':sha(z),'custody_commit':ref})
  print('custody authenticated',name,flush=True)
assert {r['bank'] for r in authenticated}==set(labels.values())
root=work/'original'
raw=load(root/'return059/RETURN059_raw_receipt.json')
mom=load(root/'return066/RETURN066_moments.json')
assert sha(root/'return059/RETURN059_raw_receipt.json')=='2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8'
assert sha(root/'return066/RETURN066_moments.json')=='b9aa229e7147b8d47d56b0cbfec290255af557cfb99260624d1c4b576369d184'
def receipt(bank):return load(banks/bank/'receipt.json')
for bank in ('091c','092a','092c','093b','099b','101b','101c','102a','102b'):assert receipt(bank)['status']=='PASS',bank
def subset(x,s,y,t):
 assert len(x)==len(y)
 if x and isinstance(x[0],int):
  assert len(x)==len(y)==2 and y[0]*s<=x[0]*t<=x[1]*t<=y[1]*s
  return 1
 return sum(subset(a,s,b,t) for a,b in zip(x,y))
records=[]
for par,sourcebank in [(0,'091c'),(1,'092a')]:
 n=518 if par==0 else 527
 arbH=load(banks/'101c'/f'H_{par}.json')
 assert arbH['parity']==par
 aux=load(banks/'092b'/'common_0.json') if par==0 else load(banks/'092c'/'odd_augmented.json')
 for suffix in ('','_alt'):
  name=f'reconstructed_matrices_{par}{suffix}.json'
  p=banks/sourcebank/name;q=banks/'092c'/name;r=banks/'093b'/name
  s,a,h=map(load,(p,q,r))
  # Each transformation must change only its declared fields.
  assert {k:v for k,v in s.items() if k not in ('T','kl','ku')}=={k:v for k,v in a.items() if k not in ('T','kl','ku')},(par,suffix,'source to auxiliary')
  assert {k:v for k,v in a.items() if k not in ('HS','Hraw')}=={k:v for k,v in h.items() if k not in ('HS','Hraw')},(par,suffix,'auxiliary to retained')
  assert a['T']==aux['scale'] and len(a['kl'])==len(a['ku'])==n
  for i in range(n):
   assert len(aux['K'][i])==n
   for j in range(n):assert aux['K'][i][j]==[a['kl'][i][j],a['ku'][i][j]]
  assert h['C']==raw['blocks'][par]['C'] and h['CS']==raw['T']
  assert F(h['target'])==(F(49,50) if par==0 else F(19,20))
  rec92=[v for v in receipt('092c')['records'] if v['input_sha256']==sha(p)]
  assert rec92 and all(v['status']=='PASS' for v in rec92)
  rr=next(v for v in receipt('093b')['records'] if v['parity']==par and v['variant']==(suffix or 'generator'))
  assert rr['source_sha256']==sha(q) and rr['sha256']==sha(r) and rr['status']=='PASS'
  assert rr['containment_transfer_verified'] and rr['positive_pivots']==87
  count=subset(arbH['Hraw'],arbH['scale'],h['Hraw'],h['HS'])
  assert count==7569
  rec={'parity':par,'variant':suffix or 'generator','source_sha256':sha(p),'auxiliary_terminal_sha256':sha(q),'retained_terminal_sha256':sha(r),'TT_U_exact_equality':True,'K_exact_equality':True,'Arb_H_containment':count}
  records.append(rec);print('literal terminal lineage PASS',par,suffix or 'generator',flush=True)
coeff=load(banks/'101b'/'receipt.json')
assert sum(subset(x,coeff['scale'],y,mom['S']) for x,y in zip(coeff['raw_coeff'],[v['raw_coeff'] for v in mom['blocks']]))==870
m=load(banks/'102b'/'moments.json');expected={(p,d['K']) for p,b in enumerate(mom['blocks']) for d in b['cutoffs']}
assert {(v['parity'],v['K']) for v in m}==expected and len(m)==len(expected)==8
n=0
for v in m:
 old=next(d for d in mom['blocks'][v['parity']]['cutoffs'] if d['K']==v['K'])
 n+=subset(v['raw_moments'],v['scale'],old['raw_moments'],mom['S'])
assert n==66816
c=load(banks/'101a'/'constants.json');d=load(banks/'101a'/'direct.json')
assert c['status']==d['status']=='PASS' and c['phase_keys']==128
old=load(root/'return086/direct_sums.json')
assert d['historical_sha256']==sha(root/'return086/direct_sums.json')
assert len(d['bands'])==2
for v in d['bands']:
 b=next(x for x in old['bands'] if x['name']==v['name'])
 assert v['count']==b['count'] and v['p']==b['p']
 # Saved Arb endpoints have been rounded outward to 1e60; historical
 # intervals are much wider, so demand containment of the saved bank.
 assert subset(v['sums'],v['scale'],b['sums'],b['scale'])==450
j=receipt('099b');e=receipt('102a')
for p in (0,1):
 r=next(v for v in e['results'] if v['parity']==p);s=next(v for v in j['results'] if v['parity']==p)
 h=F(raw['blocks'][p]['margin']);M=F(5*87,2)*F(sum(x*x for row in raw['blocks'][p]['C'] for x in row),raw['T']**2)
 assert h==F(s['retained_chart_floor']) and M==F(s['retained_mass_cap'])
 caps=[F(v['norm_upper']) for v in e['records'] if v['parity']==p]
 assert len(caps)==3 and sum(caps)==F(r['error_norm_sum'])
 rr=(F(99,100) if p==0 else F(39,40))+sum(caps)
 assert rr<1 and rr**2==F(r['true_squared_upper']) and 1-rr==F(r['relative_floor'])
 floor=(1-rr)*min(h/M,F(12,25)*F(99,1000))
 assert floor==F(r['conditional_physical_floor']) and floor>=F(3,10**59)
 print('exact physical floor binding',p,float(floor),flush=True)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'literal custody-to-terminal matrix lineage and supplemental scalar/moment/error bindings; no new mathematical promotion','archives':authenticated,'terminal_records':records,'raw_coefficients':870,'moments':66816,'direct_moments':900,'conditional_common_floor':'3/10^59'},indent=2)+'\n')
print('RETURN103-B LITERAL TERMINAL DEPENDENCY BINDING PASS',flush=True)
