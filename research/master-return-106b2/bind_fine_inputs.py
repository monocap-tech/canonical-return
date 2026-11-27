"""RETURN106-B2: resolution-correct literal finite input attachment.
Evaluate the same physical Arb balls, but export at 10^24 rather than
10^12. Compare against the authenticated original_rows at 10^16.
No widening of historical input, no overlap-based acceptance.
"""
import ast,hashlib,json,sys,subprocess,time
from pathlib import Path
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
base=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p,h):assert sha(p)==h,p
def blob(p,h):assert subprocess.check_output(['git','hash-object',str(p)],text=True).strip()==h,p
src=root/'return087/finite_audit.py'
pin(src,'cf0f1170ecaf0d5faf8b4f1560517d8f678aeba9cb000f9baf68a42fd1b0c89c')
pin(src.parent/'finite_bank.json','c6a1df49ff86806810e8763c10223b3b4b2e331b0243693802f1372eb14c443f')
p=base/'master-return-104a/regenerate_finite_arb.py'
blob(p,'802478c786d79834189fb5ed7deed811223e82a4')
blob(base/'master-return-101b/audit_raw_coefficients.py','a27db8cc14d9f389430436ec915c3150b2a18f93')
# Read the literal driver definitions, but never run main or contract.
s=src.read_text();marker="if __name__=='__main__':main()"
assert s.count(marker)==1 and s.split(marker)[1].strip()==''
legacy={'__file__':str(src),'__name__':'literal_reader_only'}
exec(compile(s.split(marker)[0],str(src),'exec'),legacy)
assert legacy['S']==10**16
old=[legacy['original_rows'](par) for par in (0,1)]
assert old[0][0]==old[1][0] and len(old[0][0])==8111
ks=old[0][0]
# Reuse the accepted physical coefficient evaluator, not its old
# output file. Change solely the final interval export resolution.
s=p.read_text();anchor='records=[];fails=[];seen=[]'
assert s.count(anchor)==1
ns={'__file__':str(p),'__name__':'physical_arb_only'}
exec(compile(s.split(anchor)[0],str(p),'exec'),ns)
assert ns['S']==10**12 and ns['Q']==10**55
assert ns['raw']==legacy['raw']
manifest=ns['manifest'];assert len(manifest)==383
# Authenticate every original source file used by the reader.
paths=['return047/RETURN047_payment_receipt.json','return063/RETURN063_constraints_receipt.json','return063/RETURN063_inverse_receipt.json','return059/RETURN059_raw_receipt.json','return062/RETURN062_bank_receipt.json','return077/proposal.json','return062/RETURN062_mixed_receipt.json','return064/RETURN064_source_receipt.json','return065/RETURN065_source_receipt.json','return067/RETURN067_source_receipt.json']+[f'return0{n}/RETURN0{n}_rows.json' for n in (68,69,71,79)]+[f'return081/rows_{i:02}.json' for i in range(24)]
for path in paths:pin(root/path,manifest[path])
assert ks==sorted(k for k,v in ns['w']['positive_even_weights'] if 116<k<=16384)
ns['S']=10**24;S=ns['S'];T=legacy['S']
records={(p,b):{'parity':p,'bank':b,'entries':0,'contained':0,'overlap':0,'disjoint':0} for p in (0,1) for b in ('early','grouped')}
fail=[];new=[[],[]];start=time.monotonic()
for i,k in enumerate(ks):
 vals=ns['mixed'](k)
 for par in (0,1):
  row=vals[par];new[par].append(row)
  assert len(row)==87
  rec=records[par,'early' if i<400 else 'grouped']
  for j,(l,u) in enumerate(row):
   a,b=map(int,old[par][1][i,j])
   assert l<=u and a<=b
   contained=a*S<=l*T<=u*T<=b*S
   overlap=a*S<=u*T and l*T<=b*S
   rec['entries']+=1;rec['contained']+=int(contained);rec['overlap']+=int(overlap);rec['disjoint']+=int(not overlap)
   if not contained and len(fail)<50:fail.append({'parity':par,'k':k,'index':j,'old':[a,b],'new':[l,u],'overlap':overlap})
 if (i+1)%512==0 or i==399:print('literal fine-scale Arb rows',i+1,'seconds',round(time.monotonic()-start,2),flush=True)
rr=list(records.values());total=sum(r['entries'] for r in rr);good=sum(r['contained'] for r in rr)
assert total==1411314
for r in rr:print('fine-scale literal containment',r,flush=True)
(out/'fine_rows.json').write_text(json.dumps({'scale':S,'ks':ks,'blocks':new},separators=(',',':'))+'\n')
receipt={'status':'PASS' if good==total else 'NOT_CERTIFIED','scope':'actual original_rows input, independent Arb physical entries at fine scale; exact integer cross multiplication','entries':total,'contained':good,'old_scale':T,'new_scale':S,'records':rr,'mismatch_samples':fail,'original_reader_sha256':sha(src),'physical_evaluator_sha256':sha(p),'original_contraction_bank_sha256':sha(src.parent/'finite_bank.json'),'fine_rows_sha256':sha(out/'fine_rows.json')}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert good==total,fail[:3]
print('RETURN106-B2 ALL LITERAL FINITE INPUT INTERVALS PASS',good,flush=True)
