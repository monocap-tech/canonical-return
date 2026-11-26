from pathlib import Path
import json, numpy as np
from fractions import Fraction as F
root=Path(__file__).resolve().parents[2]; cb=root/'ce-output/pristine/cr_int_t1d_cb'; src=root/'custody/dv/generated'; out=root/'custody/dy/generated'
# Execute the authenticated original proposal through the sampled inverse solve only.
s=(src/'bw_propose_full.py').read_text(); marker='# maximize lower bracket energy'
assert s.count(marker)==1
import os
os.chdir(cb)
exec(compile(s[:s.index(marker)],str(src/'bw_propose_full.py'),'exec'),globals())
p=np.load(src/'proposal32.npz'); assert G.shape==p['G'].shape
assert np.max(abs(G-p['G']))<1e-8
# Avoid reusing proposal variable C (cross dictionary) as the central constant.
Hinv=U.T@U*8/N;Hinv=(Hinv+Hinv.T)/2
assert np.linalg.eigvalsh(Hinv)[0]>-1e-8
baseJ=J.copy(); baseA=baseJ[d:,d:]/mu-D;baseB=K-baseJ[d:,:d]/mu
X0=np.linalg.solve(baseA,baseB); assert np.max(abs(X0-p['X']/S))<1e-8
E0=(1-eps)*np.eye(d)
def nominal(X,V):
    W=np.vstack([np.eye(d),X,-V])
    M=E0+K.T@X+X.T@K+X.T@D@X-W.T@G@W/mu
    return (M+M.T)/2
assert np.max(abs(nominal(X0,proj@np.vstack([np.eye(d),X0]))-p['central']))<1e-8
dx=json.loads((root/'custody/dx/generated/screen.json').read_text())
aq=0.0003371375047977623 # fixed head allowance recorded in DW; diagnostic only
ax=dx['physical_allowance']
rec=next(v for v in json.loads((src/'trial31/cj-summary.json').read_text()) if v['mode']==1)
delta0,t0=map(lambda x:float(F(x)),rec['parameters'][1:3])
wq,wx,we=map(lambda x:float(F(x)),rec['weights'])
assert min(delta0,t0,wq,wx,we)>0 and abs(wq+wx+we-1)<1e-12
# Proxy retains sampled inverse norm and scalar head/physical reserves.
# Omits the exact residual R, energy reserves, finite errors, and other corrections.
# In particular the proxy is NOT a certified lower bound.
rows=[]; best=None; bestarrays=None
for delta in sorted(set([delta0]+[delta0*10**v for v in np.linspace(-1.5,1.5,19)])):
    for t in [0.0,t0,2*t0]:
        gamma=delta*(1+t)
        H=G+gamma*Hinv
        pp=np.linalg.solve(H[n:,n:],H[n:,:n]); jj=H[:n,:n]-H[:n,n:]@pp;jj=(jj+jj.T)/2
        rho=(1+1/delta)*ax**2/(wx*mu); head=(1+1/delta)*aq**2/(wq*mu)
        aa=jj[d:,d:]/mu-D+rho*np.eye(m);aa=(aa+aa.T)/2
        amin=float(np.linalg.eigvalsh(aa)[0]);assert amin>0
        bb=K-jj[d:,:d]/mu
        X=np.linalg.solve(aa,bb);V=pp@np.vstack([np.eye(d),X]);W=np.vstack([np.eye(d),X,-V])
        M=nominal(X,V);inv=W.T@Hinv@W
        P=M-gamma*inv/mu-rho*X.T@X-head*np.eye(d);P=(P+P.T)/2
        ev,Q=np.linalg.eigh(P);R=Q[:,-31:];z=R[:,0]
        record=dict(delta=delta,t=t,gamma=gamma,rho=rho,head=head,proxy31=float(ev[-31]),proxy32=float(ev[-32]),nominal_on_trial=float(np.linalg.eigvalsh(R.T@M@R)[0]),nominal_rayleigh=float(z@M@z),inverse_penalty=float(gamma*(z@inv@z)/mu),physical_penalty=float(rho*np.linalg.norm(X@z)**2),hessian_min=amin)
        assert abs(record['nominal_rayleigh']-record['inverse_penalty']-record['physical_penalty']-head-record['proxy31'])<1e-8
        rows.append(record)
        if best is None or record['proxy31']>best['proxy31']:
            best=record;bestarrays=(X,V,R)
X,V,R=bestarrays
np.savez_compressed(out/'candidate.npz',X=X,V=V,Q=R,delta=best['delta'],t=best['t'])
report=dict(scope='sampled numerical proxy only; omitted exact residual and enclosure terms',gram_replay_error=float(np.max(abs(G-p['G']))),dx_proxy=dx['best']['proxy31'],head_allowance=aq,physical_allowance=ax,rows=rows,best=best)
(out/'screen.json').write_text(json.dumps(report,indent=2)+'\n')
print('DY BEST',best,flush=True)
