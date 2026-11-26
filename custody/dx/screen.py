"""DX numerical screen only: no exact lower bound or positivity claim."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import numpy as np
root=Path(__file__).resolve().parents[2]
src=root/'custody/dv/generated'
out=root/'custody/dx/generated';out.mkdir(parents=True,exist_ok=True)
files=['proposal32.npz','trial31/bs_trial.npz','merged/odd_BQ_energy.npz','merged/odd_merged_tail.npz','merged/odd_combined_projection.npz','trial31/cj-summary.json']
manifest={name.removeprefix('./'):sha for sha,name in (line.split(maxsplit=1) for line in (src/'sha256.txt').read_text().splitlines())}
for name in files:
    assert hashlib.sha256((src/name).read_bytes()).hexdigest()==manifest[name],name
p=np.load(src/files[0]); b=np.load(src/files[1]); e=np.load(src/files[2]); t=np.load(src/files[3]); y=np.load(src/files[4])
assert p['X'].shape==(352,119) and p['Q'].shape==(119,32)
d,m=119,352;n=d+m
G=p['G'];assert G.shape[0]==G.shape[1]>n and np.all(np.isfinite(G))
assert np.allclose(G,G.T,rtol=0,atol=1e-10)
mu=float(F(str(b['mu'])));eps=float(F(str(b['eps'])))
assert mu>0
K=e['K'].astype(float)/float(e['scale']);D=e['D'].astype(float)/float(e['scale'])
assert K.shape==(m,d) and D.shape==(m,m)
proj=np.linalg.solve(G[n:,n:],G[n:,:n]);J=G[:n,:n]-G[:n,n:]@proj;J=(J+J.T)/2
A=J[d:,d:]/mu-D;A=(A+A.T)/2
B=K-J[d:,:d]/mu
C=(1-eps)*np.eye(d)-J[:d,:d]/mu
def central(X):
    M=C+B.T@X+X.T@B-X.T@A@X
    return (M+M.T)/2
ae=np.linalg.eigvalsh(A);assert ae[0]>0,('nominal Hessian not positive',ae[0])
X0=np.linalg.solve(A,B);M0=central(X0)
# Compare with the original proposal and rounded coefficients before screening.
old=p['X']/float(p['scale']);oldM=p['central']
assert np.max(abs(X0-old))<1e-8
assert np.max(abs(M0-oldM))<1e-8
recs=json.loads((src/files[5]).read_text())
oldrec=next(v for v in recs if v['mode']==1)
delta=F(oldrec['parameters'][1]);wx=F(oldrec['weights'][1])
assert delta>0 and wx>0
cut=np.asarray(t['cutoff'],int);idx=np.flatnonzero(cut==16384);assert len(idx)==1
tail=F(str(t['tail'][idx[0]]));err=F(str(y['U_error']))
assert tail>=0 and err>=0
ax=float(err)+np.sqrt(float(tail))
# Fixed proxy calibrated from DV's scalar physical representation penalty.
# Omits all other variable exact reserves; do not call it a rigorous bound.
rho=(1+1/float(delta))*ax**2/(float(wx)*mu)
assert rho>0
grid=sorted(set([0.0,rho]+[rho*10**k for k in np.linspace(-3,3,25)]))
rows=[];candidates=[]
for lam in grid:
    X=np.linalg.solve(A+lam*np.eye(m),B);M=central(X)
    loss=(X-X0).T@A@(X-X0)
    assert np.max(abs(M0-M-loss))<1e-8
    proxy=M-rho*X.T@X;pev,pQ=np.linalg.eigh(proxy)
    ev,Q=np.linalg.eigh(M);Q31=Q[:,-31:];R=pQ[:,-31:]
    vals=np.linalg.eigvalsh(R.T@M@R)
    record=dict(lam=float(lam),nominal31=float(ev[-31]),nominal32=float(ev[-32]),proxy31=float(pev[-31]),nominal_on_proxy_trial=float(vals[0]),physical_norm_nominal31=float(np.linalg.norm(X@Q31,2)),physical_norm_proxy31=float(np.linalg.norm(X@R,2)),response_frobenius=float(np.linalg.norm(X)),nominal_loss_norm=float(np.linalg.norm(loss,2)))
    rows.append(record);candidates.append((X,R))
# The fixed scalar proxy has a known Loewner optimizer X_rho.
# Verify it dominates every sampled alternative numerically.
ir=grid.index(rho);Xr,_=candidates[ir]
Ar=A+rho*np.eye(m);Pr=central(Xr)-rho*Xr.T@Xr
for X,_ in candidates:
    loss=(X-Xr).T@Ar@(X-Xr)
    assert np.max(abs(Pr-(central(X)-rho*X.T@X)-loss))<1e-8
assert rows[ir]['proxy31']>=max(r['proxy31'] for r in rows)-1e-10
best=ir
X,R=candidates[best];np.savez_compressed(out/'proxy-candidate.npz',X=X,Q=R,lam=rows[best]['lam'],rho=rho,central=central(X))
report=dict(scope='numerical screen only, not an exact bound',source_commit='2e5a6cc9e09fbe69d87596c2bc86f5ee4bbe4216',mu=mu,eps=eps,hessian_min=float(ae[0]),hessian_max=float(ae[-1]),rho=rho,physical_allowance=ax,baseline_X_max_error=float(np.max(abs(X0-old))),baseline_M_max_error=float(np.max(abs(M0-oldM))),rows=rows,best=rows[best])
(out/'screen.json').write_text(json.dumps(report,indent=2)+'\n')
print('DX BASELINE',rows[0],flush=True)
print('DX BEST PROXY',rows[best],flush=True)
