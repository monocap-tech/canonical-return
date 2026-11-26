"""EL: conditional candidate32 handoff; no certified positivity.
Execute only after a positive EK numerical proxy and nominal trial.
Rebuild only 32 sampled inverse solves, not the entire 119+352 dictionary.
"""
from pathlib import Path
import hashlib,json,numpy as np
root=Path(__file__).resolve().parents[2];ek=root/'custody/ek/generated';el=root/'custody/el/generated';el.mkdir(parents=True,exist_ok=True)
for line in (ek/'sha256.txt').read_text().splitlines():
    h,name=line.split(maxsplit=1)
    with (ek/name.removeprefix('./')).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==h
report=json.loads((ek/'screen.json').read_text())
assert report['best']['proxy32']>0 and report['best']['nominal_on_trial32']>0, 'EK numerical candidate not suitable for direct promotion'
p=np.load(ek/'candidate32.npz')
X0,V0,Q0=[p[k] for k in ('X','V','Q')]
assert Q0.shape==(119,32) and X0.shape==(352,119)
assert np.max(abs(Q0.T@Q0-np.eye(32)))<1e-8
assert float(p['delta'])==report['best']['delta'] and float(p['t'])==report['best']['t']
# Reuse the authenticated numerical operator setup only. It produces no new
# screen output and performs none of the original all-dictionary inverse solves.
screen=root/'custody/ek/screen.py';src=screen.read_text()
marker='L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float);U=np.empty_like(Z)'
assert src.count(marker)==1
env=dict(__file__=str(screen),__name__='el_proposal_operator')
exec(compile(src.split(marker)[0],str(screen),'exec'),env)
S,N,d,m,xx,mask,Z,alpha,rp,C,ana,mv,LinearOperator,cg=[env[k] for k in ('S','N','d','m','xx','mask','Z','alpha','rp','C','ana','mv','LinearOperator','cg')]
assert Z.shape[1]==d+m+V0.shape[0]
assert np.all(np.isfinite(X0)) and np.all(np.isfinite(V0)) and np.all(np.isfinite(Q0))
def grid(A):
    assert np.max(abs(A))*S<2**62
    return np.rint(A*S).astype(np.int64)
q=grid(Q0);X=grid(X0);qf=q/S
v=grid(V0[:-1]@qf);g=grid(V0[-1:]@qf)
coeff=np.vstack([qf,X/S@qf,-v/S,-g/S])
rhs=Z@coeff;del Z
L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float)
UU=np.empty_like(rhs);res=[]
for i in range(32):
    y,info=cg(L,rhs[:,i],rtol=1e-9,maxiter=200);assert info==0,(i,info)
    res.append(float(np.linalg.norm(mv(y)-rhs[:,i])/max(np.linalg.norm(rhs[:,i]),1e-30)))
    UU[:,i]=y
full=np.zeros((N,32));full[mask]=UU;uc=ana(full);uu=(uc[:N//2]*np.sqrt(2)*1j).real
# Retain the existing signed-symbol cutoff 32768, independent of N.
# This is merely a trial selection: omitted inverse action needs directed bounds.
eligible=16384
idx=np.sort(np.argsort(np.sum(uu[:eligible]**2,axis=1))[-8192:])
j=2*idx+1;T=grid(uu[idx]);assert len(j)==8192 and max(j)<32768
t=grid(np.sqrt(alpha)[:,None]/(1-alpha)[:,None]*uu[rp])
assert np.array_equal(C['weights'],np.rint(alpha*1e8).astype(C['weights'].dtype))
np.savez_compressed(el/'bw_full_proposal.npz',Q=q,X=X,v0=v,g0=g,t=t,alpha=C['weights'],j=j,T=T,scale=S)
summary=dict(scope='rounded numerical proposal only; no directed certificate',dimension=32,physical_correctors=m,physical_cutoff=65536,circular_cutoff=32768,inverse_trial_modes=len(j),proxy32=report['best']['proxy32'],nominal_on_trial32=report['best']['nominal_on_trial32'],q_defect=float(np.linalg.norm(qf.T@qf-np.eye(32),2)),max_sampled_relative_residual=max(res))
(el/'proposal.json').write_text(json.dumps(summary,indent=2)+'\n')
print('EL PROPOSAL',summary,flush=True)
