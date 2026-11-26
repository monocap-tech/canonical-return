"""RETURN091-A: assemble the odd signed saddle from independently audited
tail corrections and the frozen unweighted baseline. No new numerical solve.
Shared: original baseline K, TT, U and structured retained energy H, whose
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
# The original baseline TT/U already contains full finite high band and
# old asymptotic tail. Use independently reconstructed signed corrections.
tt=sym(plus(plus(arr(base.mb['tt']),read_bank('090','restricted_bank.json',['signed_reweighting','tt'],O)),read_bank('090b','asymptotic_bank.json',['signed_reweighting','tt'],O)))
u=plus(plus(arr(base.mb['U']),read_bank('090','restricted_bank.json',['signed_reweighting','U'],O)),read_bank('090b','asymptotic_bank.json',['signed_reweighting','U'],O))
# Restoration: independent full finite high band, independent restricted
# reweighting, and direct interval for the complete strengthened tail.
c=arr(base.c['restoration_cross'])+arr(base.new['blocks'][1]['restoration_cross'])+arr(base.gap['blocks'][1]['restoration_cross'])
w=plus(plus(plus(c,read_bank('089','bank_1.json',['matrices','restoration_cross'],O)),read_bank('090','restricted_bank.json',['signed_reweighting','restoration_cross'],O)),read_bank('090b','asymptotic_bank.json',['new_restoration_cross'],O))
U=np.concatenate([u,w],axis=0)
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
receipt={'scope':'odd signed correction assembly; shared unweighted baseline and H; fixed proposals','target':str(a),'inverse_residual':str(res),'inverse_norm':str(ic),'core_radius':str(rad),'solve_reserve':str(ic*err),'margin':str(margin),'status':'PASS' if margin>0 else 'UNACCEPTED'}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
(out/'assembled.json').write_text(json.dumps({'K_scale':KS,'source_scale':O,'K':K.tolist(),'tt':tt.tolist(),'U':U.tolist(),'H_scale':S,'H':H.tolist(),'P':p.tolist()},separators=(',',':')))
assert margin>0,receipt
print('RETURN091-A ODD SIGNED CORRECTION ASSEMBLY PASS',flush=True)
