"""RETURN099-B: joint evidence binding, exact conditional physical floor.
No promotion of the analytic derivations to native Lean or a theorem.
"""
import sys,json,hashlib,math,zipfile
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();banks=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
def rec(name):return load(banks/name/'receipt.json')
with zipfile.ZipFile(root.parent/'archives/RETURN084_audit.zip') as z:
 names=[n for n in z.namelist() if n.endswith('/archive_manifest.json') or n=='archive_manifest.json'];assert len(names)==1
 data=z.read(names[0])
assert hashlib.sha256(data).hexdigest()=='38b4781ddac2b4f0437f75f8eb802e78a54c1fbc366f4690880aa73b0aa9d112'
manifest=json.loads(data)
assert len(manifest)==383
for name,h in manifest.items():assert sha(root/name)==h,name
print('original manifest files',len(manifest),flush=True)
raw=load(root/'return059/RETURN059_raw_receipt.json')
w=load(root/'return047/RETURN047_payment_receipt.json')
a=rec('093b');b=rec('094a');c=rec('095a');d=rec('095b');e=rec('096a');f=rec('096b');g=rec('097a');h=rec('097bc/early');i=rec('097bc/features')
assert all(r['status']=='PASS' for r in (a,b,c,d,e,f,g,h,i))
assert c['terminal_receipt_sha256']==sha(banks/'093b/receipt.json')
assert c['metric_receipt_sha256']==sha(banks/'094a/receipt.json')
assert b['rank']==8 and b['positive_weights']==500 and b['reference_variants']==2
assert d['negative_rows']==32 and d['finite_rows']==16302
assert {x['parity'] for x in d['records']}=={0,1}
assert all(x['target_pass'] and F(x['bound'])==F(13,25) and x['positive_pivots']==32 for x in d['records'])
assert e['chart_ranks']==[5,6] and e['retained_each']==87 and e['common_auxiliary_dimensions']==[518,519] and e['finite_complement']==16302 and e['infinite_grid_start']==32768
assert f['even_count']==16302 and f['odd_count']==16379 and not f['odd_failures'] and f['restoration_overlap']==0
assert F(f['min_even_slack'])>=0 and F(f['min_odd_slack'])>=0 and F(f['gamma_image_upper'])<F(1,1000)
assert f['restoration_receipt_sha256']==manifest['return075/RETURN075_receipt.json']
assert [(x['cutoff'],F(x['A'])) for x in f['tail']]==[(32768,F(87,500)),(65536,F(519,1000)),(262144,F(1199,1000))]
for x in f['tail']:assert F(x['symbol_lower'])>=F(x['required'])==F(1,10)+2*(F(x['A'])-F(99,1000))
# Bind every finite interval regeneration to the exact archived row labels.
expected=sorted(k for k,v in w['positive_even_weights'] if 116<k<=16384)
assert len(expected)==8111
assert g['rows']==7711 and g['entries']==g['contained']==g['overlap']==7711*174
seen=[]
for r in g['records']:
 p=root/r['file'];assert sha(p)==r['sha256']==manifest[r['file']]
 old=load(p);new=load(banks/'097a'/p.name)
 assert old['ks']==new['ks'] and old['scale']==new['scale']==10**12
 assert r['rows']==len(old['ks']) and r['entries']==r['contained']==r['overlap']==174*r['rows']
 seen+=old['ks']
 for par in (0,1):
  x=np.asarray(new['blocks'][par],dtype=object);y=np.asarray(old['blocks'][par],dtype=object)
  assert x.shape==y.shape==(r['rows'],87,2)
  assert np.all((y[:,:,0]<=x[:,:,0])&(x[:,:,0]<=x[:,:,1])&(x[:,:,1]<=y[:,:,1]))
assert seen==expected[400:]
early=load(banks/'097bc/early/early_regenerated.json')
assert early['ks']==expected[:400] and early['scale']==10**12
assert {r['parity'] for r in h['records']}=={0,1}
assert all(r['rows']==400 and r['entries']==r['contained']==r['overlap']==34800 for r in h['records'])
assert not h['mismatch_samples']
for par in (0,1):
 x=np.asarray(early['blocks'][par],dtype=object);assert x.shape==(400,87,2) and np.all(x[:,:,0]<=x[:,:,1])
features=load(banks/'097bc/features/features.json');original=load(root/'return080/features.json')
assert i['historical_features_sha256']==manifest['return080/features.json']
assert i['regenerated_features_sha256']==sha(banks/'097bc/features/features.json')
assert i['coordinates']==8192 and i['features_per_parity']==36
assert features['ks']==original['ks']==list(range(16386,32770,2))
assert features['scales']==original['scales'] and features['S']==original['S']==10**10
for par in (0,1):
 x=np.asarray(features['X'][par],dtype=object);y=np.asarray(original['X'][par],dtype=object)
 assert x.shape==y.shape==(8192,36,2)
 assert np.all((y[:,:,0]<=x[:,:,0])&(x[:,:,0]<=x[:,:,1])&(x[:,:,1]<=y[:,:,1]))
 assert i['counts'][par]['entries']==i['counts'][par]['contained']==i['counts'][par]['overlap']==294912
assert g['entries']+sum(r['entries'] for r in h['records'])+sum(r['entries'] for r in i['counts'])==2001138
print('joint finite coverage and literal grouped/feature containment PASS',2001138,flush=True)
# Independently verify the retained chart scalar floor used in the
# error envelopes via exact directed congruence of the regenerated Hraw.
# Use a deliberately crude but rigorous physical L2 comparison:
# all 87 raw basis functions have absolute value <=1 on I,
# so ||p||² <= (5/2)*87*||C/T||_F² ||x||².
results=[]
assert len(a['records'])==4 and len(c['results'])==2
for par in (0,1):
 p=next(r for r in c['results'] if r['parity']==par);target=F(49,50) if par==0 else F(19,20)
 assert F(p['modeled_squared_upper'])==target
 rcap=F(p['sqrt_cap'])+F(p['error_norm_sum'])
 assert F(p['sqrt_cap'])**2>=target and rcap**2==F(p['true_squared_upper']) and 0<rcap<1
 floor=1-rcap;assert floor==F(p['conditional_relative_form_floor'])
 C=np.asarray(raw['blocks'][par]['C'],dtype=object);T=raw['T'];assert C.shape==(87,87) and T>0
 margin=F(raw['blocks'][par]['margin']);assert margin>0
 for suffix in ('','_alt'):
  r=next(x for x in a['records'] if x['parity']==par and x['variant']==(suffix or 'generator'))
  path=banks/'093b'/f'reconstructed_matrices_{par}{suffix}.json';z=load(path)
  assert sha(path)==r['sha256'] and r['status']=='PASS' and r['containment_transfer_verified'] and r['rank_C']==r['rank_V']==87 and r['positive_pivots']==87
  assert F(z['target'])==target and z['CS']==T and z['C']==C.tolist()
  H=np.asarray(z['Hraw'],dtype=object);assert H.shape==(87,87,2) and z['HS']==10**75
  assert np.all(H[:,:,0]<=H[:,:,1]) and np.array_equal(H,H.transpose(1,0,2))
  M=C.T@(H[:,:,0]+H[:,:,1])@C;R=abs(C).T@(H[:,:,1]-H[:,:,0])@abs(C)
  m=min(int(M[j,j]-R[j,j]-sum(abs(M[j,l])+R[j,l] for l in range(87) if l!=j)) for j in range(87))
  cert=F(m,2*z['HS']*T*T)
  assert cert>=margin,(par,suffix,str(cert),str(margin))
  print('retained chart floor independently bound',par,suffix,float(cert),float(margin),flush=True)
 masscap=F(5*87,2)*F(sum(int(x)**2 for x in C.flat),T*T);assert masscap>0
 Hphysical=margin/masscap
 Dphysical=F(12,25)*F(99,1000)
 physical=floor*min(Hphysical,Dphysical)
 assert physical>0
 results.append({'parity':par,'conditional_relative_floor':str(floor),'retained_chart_floor':str(margin),'retained_mass_cap':str(masscap),'retained_physical_floor':str(Hphysical),'complement_physical_floor':str(Dphysical),'conditional_original_L2_floor':str(physical)})
 print('conditional original L2 floor',par,float(physical),flush=True)
# The chart-mass bound used sinh(5/8)<1. Prove it by a rational
# exp upper: exp(5/8)<2, from first terms and geometric tail.
x=F(5,8);n=8;s=sum(x**j/math.factorial(j) for j in range(n));assert s+x**n/math.factorial(n)/(1-x/(n+1))<2
receipt={'status':'PASS','scope':'joint exact evidence and conditional original L2 floor; analytic attachment review and permanent custody not completed; no unconditional theorem, Lean or RH','original_manifest_entries':383,'finite_interval_comparisons':2001138,'results':results,'common_conditional_original_L2_floor':str(min(F(r['conditional_original_L2_floor']) for r in results)),'receipt_hashes':{p:sha(banks/p/'receipt.json') for p in ('093b','094a','095a','095b','096a','096b','097a','097bc/early','097bc/features')}}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('RETURN099-B JOINT CONDITIONAL EVIDENCE BINDING PASS',flush=True)
