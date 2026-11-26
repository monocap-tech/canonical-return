"""RETURN091-B: assemble independent odd source baseline and signed corrections from
tail corrections and the frozen unweighted baseline. No new numerical solve.
Shared: original baseline K and structured retained energy H, whose
upstream containment checks are separate. This is not complete raw-source
regeneration, and does not establish RH.
"""
import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve(); banks=Path(sys.argv[2]).resolve(); out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
def load(p):return json.loads(p.read_bytes())
def arr(x):return np.array(x,dtype=object)
def pack(l,h):
 assert l.shape==h.shape and np.all(l<=h)
 return np.stack([l,h],axis=-1)
def iv(l,h):return pack(arr(l),arr(h))
def plus(a,b):return pack(a[:,:,0]+b[:,:,0],a[:,:,1]+b[:,:,1])
def neg(a):return pack(-a[:,:,1],-a[:,:,0])
def diff(a,b):return plus(a,neg(b))
def scale(a,fac):
 fac=F(fac)
 if fac<0:return scale(neg(a),-fac)
 p,q=fac.numerator,fac.denominator
 return pack(a[:,:,0]*p//q,-((-a[:,:,1]*p)//q))
def sym(a):
 assert a.shape[0]==a.shape[1]
 return pack(np.maximum(a[:,:,0],a[:,:,0].T),np.minimum(a[:,:,1],a[:,:,1].T))
def tr(a):return a.transpose(1,0,2)
def right(a,x):
 p=np.maximum(x,0);n=np.minimum(x,0);l,h=a[:,:,0],a[:,:,1]
 return pack(l@p+h@n,h@p+l@n)
def read_bank(folder,file,key,scaleout):
 d=load(banks/folder/file); a=d
 for k in key:a=a[k]
 return scale(arr(a),F(scaleout,d['scale']))
def enc(a,s):return math.floor(a[0]*s),math.ceil(a[1]*s)
# Archived setup is a source of *baseline* inputs, not terminal arithmetic.
sys.path.insert(0,str(root/'return083'))
sys.argv=[str(root/'return083/upgrade.py')]
import setup as base
assert not base.alt
O=10**12; KS=10**8; XS=10**7; YS=10**12; VS=10**12
assert base.O==O
lam=F(174,1199);theta=1-lam;theta2=F(115,173)
# Original augmented saddle with signed inverse metric and restoration ridge.
l=np.zeros((527,527),object);h=l.copy()
K0=scale(iv(base.KL,base.KU),F(KS,base.T))
l[:519,:519]=K0[:,:,0];h[:519,:519]=K0[:,:,1]
W=scale(arr(base.W),F(KS,base.WS))
l[:519,519:]=W[:,:,0];h[:519,519:]=W[:,:,1]
l[519:,:519]=W[:,:,0].T;h[519:,:519]=W[:,:,1].T
for i in range(8):
 for j in range(8):
  a,b=map(F,base.r['gram'][i][j]);v=F(base.r['R'][i][j])
  l[519+i,519+j],h[519+i,519+j]=enc((a+v,b+v),KS)
Kbase=pack(l,h)
# Independent tail Gram corrections: retain a common scale before subtraction.
g0=read_bank('088','gram_262144.json',['gram'],10**40)
g2=read_bank('088','gram_65536.json',['gram'],10**40)
K=sym(plus(Kbase,scale(plus(scale(g0,-theta),scale(diff(g2,g0),-theta2)),F(KS,10**40))))
# RETURN087 banks are read directly from the authenticated archive, not
# regenerated or substituted by parent values. Bind the explicit archived parity and row-range sections, then require
# full three-matrix containment against their parent sections.
r87=root/'return087'
pins={'finite_bank.json':'c6a1df49ff86806810e8763c10223b3b4b2e331b0243693802f1372eb14c443f',
      'band_0.json':'624d811d17b570bdc58d39d7e47350bad0f91524443e8e2fad0cba91a7c6db58',
      'band_1.json':'89805afe0b4fdcef6a47847490d6940b75e3d95b0f5203cdc74f087098c319c1'}
for name,sha in pins.items():
 assert hashlib.sha256((r87/name).read_bytes()).hexdigest()==sha,name
keys=('tt','U','restoration_cross')
shapes={'tt':(87,87,2),'U':(519,87,2),'restoration_cross':(8,87,2)}
def contained(a,b):
 return np.all(a[:,:,0]>=b[:,:,0]) and np.all(a[:,:,1]<=b[:,:,1])
failures=[]
def checked(bank,sc,parent,label):
 assert all(k in bank and k in parent for k in keys),(label,list(bank),list(parent))
 result={}
 for k in keys:
  a=scale(arr(bank[k]),F(O,sc));p=arr(parent[k])
  assert a.shape==p.shape==shapes[k],(label,k,a.shape,p.shape)
  inside=(a[:,:,0]>=p[:,:,0])&(a[:,:,1]<=p[:,:,1])
  overlap=(a[:,:,0]<=p[:,:,1])&(p[:,:,0]<=a[:,:,1])
  print('RETURN087 containment',label,k,int(inside.sum()),inside.size,'overlap',int(overlap.sum()),'sample',a[0,0].tolist(),p[0,0].tolist(),flush=True)
  failures.append((label,k,int(inside.sum()),inside.size))
  result[k]=a
 return result
f87=load(r87/'finite_bank.json');b87=load(r87/'band_1.json')
# Bind using the audited bank's parity and ordered row ranges, not a
# search for a favorable numerical containment match.
odd=[(i,v) for i,v in enumerate(f87['banks']) if v['parity']==1]
assert len(odd)==3
print('RETURN087 odd partitions',[(i,v['start'],v['end']) for i,v in odd],flush=True)
assert [(v['start'],v['end']) for _,v in odd]==[(0,1680),(1680,2192),(2192,8111)]
ie,eb=odd[1];il,lb=odd[2]
assert ie==4 and il==5
ep=base.c['blocks'][1]
print('RETURN079 finite keys',list(ep['finite']),'scales',ep['scales'],'root scale',base.c['scale'],flush=True)
early_parent=dict(ep['finite'],restoration_cross=base.c['restoration_cross'])
assert base.c['scale']==O
early=checked(eb['cumulative'],f87['scale'],early_parent,'cumulative 2192 actual rows')
# Compare against original serialized input, not mutable objects imported
# from a legacy assembly module. Print archived verifier's binding code
# and exact metadata before evaluating all three banks.
for f in (root/'return083/setup.py',r87/'verify.py'):
 print('authenticated source binding',str(f.relative_to(root)),flush=True)
 for line in f.read_text().splitlines():
  if any(t in line for t in ('contract.json','contract_alt.json','band.json','finite_bank','increment','new=','gap=','mb=','scale','base.new')):print(line,flush=True)
raw_late=load(root/'return081/contract.json')
raw_band=load(root/'return080/band.json')
assert raw_late['scale']==raw_band['scale']==O
assert len(raw_late['ks'])==5919 and raw_late['ks'][0]==4544 and raw_late['ks'][-1]==16384
assert len(raw_band['ks'])==8192 and raw_band['ks'][0]==16386 and raw_band['ks'][-1]==32768
for label,raw,loaded in [('late',raw_late,base.new),('band',raw_band,base.gap)]:
 print('raw/legacy equality',label,{k:bool(np.array_equal(arr(raw['blocks'][1][k]),arr(loaded['blocks'][1][k]))) for k in keys},flush=True)
late=checked(lb['increment'],f87['scale'],raw_late['blocks'][1],'increment 5919 actual rows')
assert b87['record']['parity']==1
band=checked(b87['bank'],O,raw_band['blocks'][1],'8192 weighted feature rows')
(out/'containment.json').write_text(json.dumps(failures,indent=2)+'\\n')
assert all(n==d for _,_,n,d in failures),failures
pe=['banks',str(ie),'cumulative'];pl=['banks',str(il),'increment'];pb=['bank']
# The independently audited source partition is 8111 actual rows,
# 8192 weighted-model rows, full finite high band, and infinite tail.
# The restricted band is a signed correction to the finite high band only.
source={}
for k in keys:
 full=read_bank('089','bank_1.json',['matrices',k],O)
 restricted=read_bank('090','restricted_bank.json',['signed_reweighting',k],O)
 tail=(read_bank('090b','asymptotic_bank.json',['new_restoration_cross'],O) if k=='restoration_cross'
       else scale(read_bank('090b','asymptotic_bank.json',['old',k],O),lam))
 source[k]=plus(plus(plus(plus(plus(early[k],late[k]),band[k]),full),restricted),tail)
tt=sym(source['tt']);u=source['U'];w=source['restoration_cross']
U=np.concatenate([u,w],axis=0)
# Compare to the previously accepted 091-A assembled source; require overlap
# entry by entry but do not intersect or reuse its values in the new solve.
prior=load(banks/'091a'/'assembled.json');assert prior['source_scale']==O
for k,m in [('tt',tt),('U',U)]:
 p=arr(prior[k]);assert p.shape==m.shape
 overlap=(m[:,:,0]<=p[:,:,1])&(p[:,:,0]<=m[:,:,1])
 print('independent baseline source overlap',k,int(overlap.sum()),overlap.size,flush=True)
 assert np.all(overlap),k
(out/'source_partition.json').write_text(json.dumps({'source_rows':[2192,5919,8192],'paths':{'early':pe,'late':pl,'band':pb},'pins':pins,'even_bank_authenticated_not_assembled':True},indent=2)+'\n')
assert tt.shape==(87,87,2) and U.shape==(527,87,2) and K.shape==(527,527,2)
assert np.all(tt[:,:,0]<=tt[:,:,1]) and np.all(U[:,:,0]<=U[:,:,1])
print('assembled signed odd saddle',K.shape,U.shape,tt.shape,flush=True)
# Fixed old inverse, solve and congruence proposals; check all with exact
# Python integer products and directed endpoints, no floating optimization.
ref=load(root/'return083/upgrade.json')
X=arr(ref['X']);Y=arr(ref['Y']);V=arr(ref['V'])
assert X.shape==(527,527) and Y.shape==(527,87) and V.shape==(87,87)
assert np.array_equal(X,X.T)
kx=right(K,X);I=KS*XS*np.eye(527,dtype=object)
res=F(int(max(np.sum(np.maximum(abs(I-kx[:,:,0]),abs(I-kx[:,:,1])),axis=1))),KS*XS)
assert res<1,('inverse residual',str(res))
ic=F(int(max(np.sum(abs(X),axis=1))),XS)/(1-res)
print('inverse residual',float(res),'inverse norm upper',float(ic),flush=True)
ky=right(K,Y);uy=right(tr(U),Y);yky=right(tr(ky),Y)
D=O*KS*YS**2
core=sym(diff(plus(scale(tt,D//O),scale(yky,D//(KS*YS**2))),scale(plus(uy,tr(uy)),D//(O*YS))))
E=diff(scale(U,KS*YS),scale(ky,O))
err=F(sum(int(z)**2 for z in np.maximum(abs(E[:,:,0]),abs(E[:,:,1])).flat),(O*KS*YS)**2)
m=(core[:,:,0]+core[:,:,1])//2
rad=F(int(max(np.sum(np.maximum(m-core[:,:,0],core[:,:,1]-m),axis=1))),D)
S=base.S
# Symmetric inverse need not be positive: operator norm controls signed error.
extra=math.ceil(F(4,5)*(rad+ic*err)*S)
p=pack(4*m*S//(5*D),-((-4*m*S)//(5*D)))
p+=np.stack([extra*np.eye(87,dtype=object)]*2,axis=-1)
# Reuse the archived structured-energy lower construction only.
s=(root/'return081/diagnostic.py').read_text()
env={'base':base.pay.base,'np':np,'math':math,'S':S}
exec(s[s.index('def structured_H'):s.index('def dot(')],env)
hl,hu,_,_=env['structured_H'](1); H=iv(hl,hu)
a=F(19,20)
# Integer numerator of aH-P at scale S*a.denominator.
C=diff(scale(H,a.numerator),scale(p,a.denominator))
C=sym(right(tr(right(C,V)),V))
lower=C[:,:,0];upper=C[:,:,1]
margin=min(int(lower[i,i])-sum(max(abs(int(lower[i,j])),abs(int(upper[i,j]))) for j in range(87) if i!=j) for i in range(87))
margin=F(margin,S*a.denominator*VS**2)
print('target',float(a),'positive congruence margin',float(margin),'radius',float(rad),'solve reserve',float(ic*err),flush=True)
receipt={'scope':'odd independently reconstructed source baseline; shared original auxiliary K and H; fixed proposals','target':str(a),'inverse_residual':str(res),'inverse_norm':str(ic),'core_radius':str(rad),'solve_reserve':str(ic*err),'margin':str(margin),'status':'PASS' if margin>0 else 'UNACCEPTED'}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
(out/'assembled.json').write_text(json.dumps({'K_scale':KS,'source_scale':O,'K':K.tolist(),'tt':tt.tolist(),'U':U.tolist(),'H_scale':S,'H':H.tolist(),'P':p.tolist()},separators=(',',':')))
assert margin>0,receipt
print('RETURN091-B ODD INDEPENDENT SOURCE BASELINE PASS',flush=True)
