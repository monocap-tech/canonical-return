from pathlib import Path
import os,json,numpy as np
root=Path(__file__).resolve().parents[2];cb=root/'ce-output/pristine/cr_int_t1d_cb';src=root/'custody/dv/generated';out=root/'custody/dz/generated'
s=(src/'bw_propose_full.py').read_text()
marker='# maximize lower bracket energy';tail='# numerical joint inverse of actual rounded proposed physical representative'
assert s.count(marker)==s.count(tail)==1
os.chdir(cb)
exec(compile(s[:s.index(marker)],str(src/'bw_propose_full.py'),'exec'),globals())
frozen=np.load(src/'proposal32.npz')
assert np.max(abs(G-frozen['G']))<1e-8
candidate=np.load(root/'custody/dy/generated/candidate.npz')
report=json.loads((root/'custody/dy/generated/screen.json').read_text())
X=candidate['X'];Vfull=candidate['V'];Qtrial=candidate['Q']
assert X.shape==(352,119) and Vfull.shape==(G.shape[0]-n,119) and Qtrial.shape==(119,31)
assert np.max(abs(Qtrial.T@Qtrial-np.eye(31)))<1e-8
assert float(candidate['delta'])==report['best']['delta'] and float(candidate['t'])==report['best']['t']
W=np.vstack([np.eye(d),X,-Vfull])
M=(1-eps)*np.eye(d)+K.T@X+X.T@K+X.T@D@X-W.T@G@W/mu;M=(M+M.T)/2
target=float(np.linalg.eigvalsh(Qtrial.T@M@Qtrial)[0])
assert abs(target-report['best']['nominal_on_trial'])<1e-8
ntrial=31;q=np.rint(Qtrial[:,::-1]*S).astype(np.int64);qf=q/S
Xgrid=np.rint(X*S).astype(np.int64);V=Vfull@qf
vgrid=np.rint(V[:-1]*S).astype(np.int64);ggrid=np.rint(V[-1:]*S).astype(np.int64)
exec(compile(s[s.index(tail):],str(src/'bw_propose_full.py'),'exec'),globals())
p=np.load('bw_full_proposal.npz');assert np.array_equal(p['Q'],q) and np.array_equal(p['X'],Xgrid)
record=dict(scope='numerical proposal only',proxy=report['best']['proxy31'],nominal_on_trial=target,rounded_q_defect=float(np.linalg.norm(qf.T@qf-np.eye(31),2)),sampled_parameter_t=report['best']['t'])
(out/'proposal-check.json').write_text(json.dumps(record,indent=2)+'\n')
print('DZ PROPOSAL',record,flush=True)
