from pathlib import Path
import json,hashlib
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
exec((ROOT.parent/'cr_alp_054/check_saved.py').read_text().split('result={}')[0].replace('ROOT=Path(__file__).resolve().parent',''))
import sys
suffix=sys.argv[1] if len(sys.argv)>1 else ''
result={}
for k in (1,2):
 raw=(ROOT/f'correction_{k}.json').read_bytes();c=json.loads(raw);r=json.loads((ROOT/f'result{suffix}_{k}.json').read_text());assert hashlib.sha256(raw).hexdigest()==r['correction_sha256']; assert hashlib.sha256((ROOT.parent/'cr_alp_055'/f'bank_{k}.json').read_bytes()).hexdigest()==c['bank_sha256']
 S=int(r['scale']);K,G,H=map(load,[r['K'],r['Gin'],r['Gout']]);s=S//int(c['scale']);assert K==[[(l*s,h*s) for l,h in row] for row in load(c['K'])]; assert G==[[(l*s,h*s) for l,h in row] for row in load(c['G'])];assert ld(K,F(1,16),S)==44-k;assert upper(G,S)<F(1,40000)
 B=[[(floor(F(200,51)*(l+a)),ceil(F(200,51)*(h+b))) for (l,h),(a,b) in zip(row,rr)] for row,rr in zip(G,H)];C=sub(K,B);nw=quad(H,r['norm_witness']['z'],S);cw=quad(C,r['comparison_witness']['z'],S);assert nw[0]>F(1,100)
 sh=None
 for q in map(F,['1/100','1/1000','1/10000','0']):
  try:ld(C,q,S);sh=q;break
  except AssertionError:pass
 status='PASS' if sh is not None else 'FALSIFIED' if cw[1]<0 else 'INCONCLUSIVE'
 print(k,r['end'],'norm lower',float(nw[0]),'comparison upper',float(cw[1]),'scalar',status,'shift',sh,flush=True)
 result[str(k)]={'cutoff':r['end'],'norm_lower':str(nw[0]),'comparison_upper':str(cw[1]),'scalar_status':status,'shift':str(sh),'corrected_margin':'1/16','inband_gram_upper':str(upper(G,S))}
(ROOT/f'verification{suffix}.json').write_text(json.dumps(result,indent=2)+'\n')
