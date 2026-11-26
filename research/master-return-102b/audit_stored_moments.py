"""RETURN102-B: strict physical raw moment containment at historical anchors.
Reuse the accepted independent Arb coefficient evaluator, not its
serialized rounded intervals. Authenticate and compare actual source moments.
"""
import sys,json,hashlib
from pathlib import Path
# The accepted evaluator defines coefficient(), ids, mom and enc before
# its final driver. Execute only those definitions and input authentication.
p=Path(__file__).resolve().parents[1]/'master-return-101b/audit_raw_coefficients.py'
src=p.read_text();anchor='records=[];allrows=[];fail=[]'
assert src.count(anchor)==1
env={'__file__':str(p),'__name__':'independent_coefficients'}
exec(compile(src.split(anchor)[0],str(p),'exec'),env)
from flint import arb,ctx
from fractions import Fraction as F
ctx.prec=384
root=env['root'];out=env['out'];raw=env['raw'];mom=env['mom'];ks0=env['ks0'];coefficient=env['coefficient'];enc=env['enc'];A=env['A'];pi=env['pi']
assert mom['S']==10**55 and mom['T']==raw['T'] and mom['N']==32
Q=mom['S']
lead=env['load']('return060/RETURN060_leading_receipt.json','e9001e8dfad3e1a0fdedbaba9d27c258afb241134d16f98bbff7e6477007fe57')
records=[];fail=[];banks=[]
def check(v,old,tag):
 l,u=enc(v);a,b=old;assert a<=b and l<=u
 if not a<=l<=u<=b:
  if len(fail)<20:fail.append([*tag,[l,u],old])
  return 0
 return 1
for par in (0,1):
 ids=([(0,0)] if par==0 else [])+[(0,k) for k in ks0]+([(1,0)] if par else [])+[(0,k) for k in raw['ks']]
 assert len(ids)==87
 cf=[coefficient(par,*x) for x in ids]
 bb=lead['blocks'][par]
 assert len(bb['raw_boundary'])==87
 count=sum(check(row[1],old,['boundary',par,i]) for i,(row,old) in enumerate(zip(cf,bb['raw_boundary'])))
 records.append({'parity':par,'kind':'leading_boundary','entries':87,'contained':count})
 # The historical odd constant is 2c+d-4b (RETURN066 replay).
 if par:
  assert len(bb['raw_constant'])==87
  count=sum(check(2*row[3]+row[4]-4*row[1],old,['odd_constant',i]) for i,(row,old) in enumerate(zip(cf,bb['raw_constant'])))
  records.append({'parity':par,'kind':'leading_constant','entries':87,'contained':count})
 cutoffs=mom['blocks'][par]['cutoffs'];assert len(cutoffs)>0
 for dat in cutoffs:
  K=dat['K'];R=2*pi*K/5;pr=-A(F(1,4))/R**2
  # The exact integer trigonometric ratio avoids cancellation of pi.
  rat=[A(F(k*k,K*K)) if not typ else pr for typ,k in ids]
  old=dat['raw_moments'];assert len(old)==32
  count=0;result=[]
  for n,ref in enumerate(old):
   assert len(ref)==3 and all(len(x)==87 for x in ref)
   v=[[row[j]*t**n+(row[4]*pr**n/2 if j==3 else arb(0)) for row,t in zip(cf,rat)] for j in (1,2,3)]
   for j in range(3):
    for i,x in enumerate(v[j]):count+=check(x,ref[j][i],['moment',par,K,n,j,i])
   result.append([[enc(x)for x in row]for row in v])
  rec={'parity':par,'kind':'raw_moments','K':K,'entries':32*3*87,'contained':count}
  records.append(rec);banks.append({'parity':par,'K':K,'scale':Q,'raw_moments':result})
  print('literal stored moment containment',rec,flush=True)
ok=all(r['contained']==r['entries'] for r in records)
(out/'receipt.json').write_text(json.dumps({'status':'PASS' if ok else 'NOT_CERTIFIED','scope':'independent physical raw moments and leading endpoint coefficient containment, not error norm or terminal recomputation','records':records,'mismatches':fail,'coefficient_source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2)+'\n')
(out/'moments.json').write_text(json.dumps(banks,separators=(',',':'))+'\n')
assert ok,fail[:3]
print('RETURN102-B LITERAL STORED MOMENTS PASS',sum(r['entries']for r in records),flush=True)
