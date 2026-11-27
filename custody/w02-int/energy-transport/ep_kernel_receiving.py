from pathlib import Path
import numpy as np,json,zipfile,io,hashlib
from fractions import Fraction as F
root=Path(__file__).resolve().parent
pins={'EO':'7e519ccdca9440b65103fffa20b24d8e6edeba9b','EN':'d84732c72265bfa30f86463f8108c9e39f62880e','EG':'06805c7f3bb20286f42ef00627937990cb9bf308','DZ':'036d47d3710631055e0cc6e7458af7ab339ae688'}
def load(s,n):
 with zipfile.ZipFile(root/f'EP-parent-{s}-evidence.zip') as z:
  m=z.read('sha256.txt');assert hashlib.sha1(b'blob '+str(len(m)).encode()+b'\0'+m).hexdigest()==pins[s]; hashes={p.removeprefix('./'):d for d,p in (line.split(maxsplit=1) for line in m.decode().splitlines())};raw=z.read(n);assert hashlib.sha256(raw).hexdigest()==hashes[n];return dict(np.load(io.BytesIO(raw),allow_pickle=False))
ka=np.load(root/'EP-receiver-kernel-proposal.npz',allow_pickle=False);ids=ka['indices'];T=ka['coefficient'].astype(float)/float(ka['coefficient_scale']);Yk=ka['action_center'];h=load('EG','odd_head_full.npz');assert np.array_equal(ka['action_indices'],h['indices'])
phys=load('EN','merged/odd_merged_corrector.npz');U=phys['coefficient'].astype(float)/float(phys['coefficient_scale']);e=load('EN','merged/odd_BQ_energy.npz');p=load('EN','merged/odd_combined_projection.npz');K=e['K'].astype(float)/float(e['scale']);D=e['D'].astype(float)/float(e['scale']);Y=p['Y']/float(p['Yscale'])
rr=lambda ids:(.25+(np.pi/4)**2)/(.25+(ids*np.pi/4)**2)
def pairing(ids,T,Y):return np.sqrt(2)*T.T@(Y[(ids-1)//2]-rr(ids)[:,None]*Y[:1])
oldE=pairing(phys['indices'],U,Y);ctrlK=float(np.max(abs(oldE[:,:119]-K)));ctrlD=float(np.max(abs(oldE[:,119:]-D)));assert max(ctrlK,ctrlD)<1e-9
cross=pairing(ids,T,Y);self0=pairing(ids,T,Yk);sym=float(np.max(abs(self0-self0.T)));selfE=(self0+self0.T)/2;assert sym<1e-8
cross_check=float(np.max(abs(pairing(phys['indices'],U,Yk)-cross[:,119:].T)))
a=load('EO','candidate32.npz');ap=load('EO','candidate32_positive_t.npz');q=np.column_stack([a['Q'][:,-1],ap['Q'][:,-1]]);x=np.column_stack([a['X']@q[:,0],ap['X']@q[:,1]]);W=np.vstack([q,x]);Ec=cross@W
out={'scope':'Numerical kernel mixed energy and action centers; no new directed interval bounds or inverse metric comparison.','control_K_max_abs':ctrlK,'control_D_max_abs':ctrlD,'self_symmetry_defect':sym,'mixed_reciprocity_max_abs':cross_check,'self_energy_eigenvalues':np.linalg.eigvalsh(selfE).tolist(),'self_energy_norm':float(np.linalg.norm(selfE,2)),'head_cross_norm':float(np.linalg.norm(cross[:,:119])),'physical_cross_norm':float(np.linalg.norm(cross[:,119:])),'frozen_ray_energy_coupling_norms':np.linalg.norm(Ec,axis=0).tolist()}
(root/'EP-kernel-energy.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(root/'EP-kernel-energy.npz',K=cross[:,:119],Dcross=cross[:,119:],Dnew=selfE,ray_coupling=Ec);print(json.dumps(out,indent=2))
# Independently evaluate the frozen EO exceptional representatives and kernel response.
pins.update({'DV':'dc58e2c0003c4d47ccaec7ab2251a6af4890f48b','CR':'85b77e4beaddbeba6c69a40087db283bc52c182f','DA':'5bd5ba86d8454e6c730c03f4ed890ef3587957f1','CK':'eb81cbe578e5efb3745daf6229ac2e089b3704ff','CQ':'77e08780d29f3ffe9cc6c0e8785598ab5a94557a','DG':'40392aeea044e5d05e81495dc0ccbc58249f6b5e','DK':'78caa974fe4df2e8329e1a6719861e87f9445eea'})
V=[];C=[];js=weights=None
for stage in ('DV','DZ','CR','DA','CK','CQ','DG','DK'):
 with zipfile.ZipFile(root/f'EP-parent-{stage}-evidence.zip') as z: names=sorted(n for n in z.namelist() if n.endswith('/bs_trial.npz'))
 for n in names:
  t0=load(stage,n)
  if js is None:js=t0['js'];weights=t0['a']
  assert np.array_equal(js,t0['js']) and np.array_equal(weights,t0['a']);V.append(t0['v0']/1e12);C.append(-t0['c'][:-len(weights)]/1e12)
V=np.column_stack(V);C=np.column_stack(C);u,s,ht=np.linalg.svd(V,full_matrices=False);rj=js[-len(weights):];jf=js[:-len(weights)];alpha=weights/1e8
src=(root/'ep_expanded_exceptional_replay.py').read_text();from scipy.sparse.linalg import LinearOperator,cg
exec(compile(src[src.index('class Op:'):src.index('op=Op(32768)')],str(root/'ep_expanded_exceptional_replay.py'),'exec'))
def report(stage):
 with zipfile.ZipFile(root/f'EP-parent-{stage}-evidence.zip') as z:
  m=z.read('sha256.txt');hs={p.removeprefix('./'):d for d,p in (l.split(maxsplit=1) for l in m.decode().splitlines())};assert hashlib.sha1(b'blob '+str(len(m)).encode()+b'\0'+m).hexdigest()==pins[stage];raw=z.read('extension.json');assert hashlib.sha256(raw).hexdigest()==hs['extension.json'];return json.loads(raw)
eg=report('EG');en=report('EN');old=[p['Y'][:,119:471]];new=[]
for b0 in eg['new_chunks']:old.append(load('EG',f"odd_U_{b0['lo']}_{b0['hi']}.npz")['Y'])
for b0 in en['new_chunks']:new.append(load('EN',f"new/odd_U_{b0['lo']}_{b0['hi']}.npz")['Y'])
Ya=np.column_stack([h['Y'],np.vstack(old),np.vstack(new)])/1e12;assert np.array_equal(Ya[:8192],Y)
# Fixed EO exceptional coefficients, rank-115 dictionary reconstruction; no exceptional reoptimization.
ve=np.column_stack([a['V']@q[:,0],ap['V']@q[:,1]]);rank=115;ec=((C+.5e-12)@ht[:rank].T/s[:rank])@(u[:,:rank].T@ve[:-1]);op=Op(131072);Z=op.neu(h['indices'],np.column_stack([Ya@W,Yk]));Z[:,:2]-=op.sin(jf,ec)+np.sinh(op.xx/2)[:,None]*ve[-1:];sol=[];res=[]
for col in Z.T:
 y0,r0=op.solve(col);sol.append(y0);res.append(r0)
I=np.column_stack(sol);G=Z.T@I*op.dx;G=(G+G.T)/2;H2=I.T@I*op.dx
b0=load('DZ','trial31/bs_trial.npz');mu=float(F(str(b0['mu'])));eps=float(F(str(b0['eps'])));weights0=list(map(F,'592921980/1249999999 2077677507/9999999992 315059687/9999999992 109269713/4999999996 20840557/1249999999 277829459/9999999992 114897809/9999999992 266801263/9999999992 237898775/9999999992 210888273/9999999992 33960425/4999999996 55140225/4999999996 111464579/9999999992 13701323/2499999998 116271951/9999999992 137267045/9999999992 340209197/9999999992 153951925/9999999992 4762100/1249999999 5920425/2499999998 54089427/2499999998'.split()));aq=F('2235124560084063602688731000016297644723/500000000000000000000000000000000000000000000');er=F('27795457820755002025071/125000000000000000000000000000');assert sum(weights0)==1
E0=(1-eps)*q.T@q+x.T@K@q+q.T@K.T@x+x.T@D@x;E=np.block([[E0,Ec.T],[Ec,selfE]]);records=[]
for ii,aa in enumerate((a,ap)):
 delta=float(aa['delta']);t=float(aa['t']);mass=float(aa['new_block_weight']);ws=[(1-mass)*float(w) for w in weights0[:-1]]+[mass,(1-mass)*float(weights0[-1])];rho=np.empty(368)
 for j,b in enumerate(en['physical_blocks']):rho[b['lo']:b['hi']]=(1+1/delta)*float(F(b['allowance']))**2/ws[j+1]/mu
 head=(1+1/delta)*(float(aq)**2/ws[0]+float(er)**2/ws[-1])/mu;P=E-(G+delta*(1+t)*H2)/mu;P[:2,:2]-=x.T@(rho[:,None]*x)+head*q.T@q;P=(P+P.T)/2;A=-P[2:,2:];bb=P[2:,ii];ev=np.linalg.eigvalsh(A);base=P[ii,ii];rec={'variant':ii,'baseline':float(base),'hessian_min':float(ev[0]),'hessian_max':float(ev[-1]),'coupling_norm':float(np.linalg.norm(bb)),'reserves':[]}
 for reserve in (0.,.01,.1,1.,10.,100.):
  if ev[0]+reserve<=0:continue
  yy=np.linalg.solve(A+reserve*np.eye(16),bb);gain=float(bb@yy);rec['reserves'].append({'artificial_reserve':reserve,'gain':gain,'optimized':float(base+gain),'coefficient_norm':float(np.linalg.norm(yy))})
 records.append(rec)
 np.savez_compressed(root/f'EP-kernel-receiving-block-{ii}.npz',A=A,b=bb,baseline=base)
out.update(scope='Numerical kernel receiving with frozen EO head, physical and exceptional coefficients; new kernel coordinates only optimized. Inherited penalties retained; new directed error/tail penalties absent. Artificial reserves exploratory.',sample_count=op.N,solve_count=len(res),max_solve_residual=max(res),rank=rank,records=records);(root/'EP-kernel-receiving.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
