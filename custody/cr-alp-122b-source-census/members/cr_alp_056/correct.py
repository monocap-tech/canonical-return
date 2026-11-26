from core import *
for k in (1,2):
 raw=(ROOT.parent/'cr_alp_055'/f'bank_{k}.json').read_bytes(); assert hashlib.sha256(raw).hexdigest()=={1:'1ec38f55beb2d59e1467be94908fe1ca33d45d24e7cfab310084775e65f34bdb',2:'47a10b15f4bdd68e5b84fb98303455df3bd3b638eeed06c76b2f88bfc184a373'}[k]
 b=json.loads(raw); P,RF,K=[load(b[a],int(b['scale'])) for a in ('P','R','K')]; gamma=F(72,25);delta=F(21,8);x=gamma/delta; D=scale(P,-1/delta); Z0=zero(RF[0].shape);Z1=scale(RF,1/gamma);t0,t1=F(1),x
 for m in range(2,13):
  t2=2*x*t1-t0; Z2=sub(scale(add(mul(D,Z1),scale(RF,1/delta)),2*t1/t2),scale(Z0,t0/t2));Z0,Z1,t0,t1=Z1,Z2,t1,t2
 L=add(P,scale(identity(len(P[0])),gamma));R=sub(RF,mul(L,Z1));cross=mul(trans(RF),Z1);Kp=symmetric(add(sub(sub(K,cross),trans(cross)),mul(trans(Z1),mul(L,Z1))));G=symmetric(mul(trans(R),R));assert upper(G)<F(1,40000);assert cert(Kp,F(1,16))==len(K[0]);
 print(k,'K margin 1/16; trial residual Gram upper',float(upper(G)),flush=True)
 (ROOT/f'correction_{k}.json').write_text(json.dumps({'scale':str(S),'Z':bounds(Z1),'K':bounds(Kp),'G':bounds(G),'bank_sha256':hashlib.sha256(raw).hexdigest()})+'\n')
