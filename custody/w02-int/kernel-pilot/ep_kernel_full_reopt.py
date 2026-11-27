from pathlib import Path
import numpy as np,json,time
root=Path(__file__).resolve().parent
src=(root/'ep_kernel_receiving.py').read_text()
exec(compile(src[:src.index("a=load('EO','candidate32.npz')")], 'energy input', 'exec'))
exec(compile(src[src.index('pins.update('):src.index('# Fixed EO exceptional coefficients')], 'operator input', 'exec'))
a=load('EO','candidate32.npz');ap=load('EO','candidate32_positive_t.npz');q=np.column_stack([a['Q'],ap['Q']]);x=np.column_stack([a['X']@a['Q'],ap['X']@ap['Q']]);W=np.vstack([q,x]);ve=np.column_stack([a['V']@a['Q'],ap['V']@ap['Q']]);rank=115
op=Op(131072);basis=(C+.5e-12)@ht[:rank].T/s[:rank];ec=basis@(u[:,:rank].T@ve[:-1]); Z=np.empty((len(op.xx),64+16+368+rank+1));Z[:,:64]=op.neu(h['indices'],Ya@W)-op.sin(jf,ec)-np.sinh(op.xx/2)[:,None]*ve[-1:];Z[:,64:80]=op.neu(h['indices'],Yk)
for j in range(0,368,16):Z[:,80+j:80+min(j+16,368)]=op.neu(h['indices'],Ya[:,119+j:119+min(j+16,368)])
Z[:,448:448+rank]=op.sin(jf,basis);Z[:,-1]=np.sinh(op.xx/2);del Ya
I=np.empty_like(Z);res=[];start=time.time()
for i in range(Z.shape[1]):
 I[:,i],r0=op.solve(Z[:,i]);res.append(r0)
 if i%16==15 or i==Z.shape[1]-1:print('solved',i+1,'of',Z.shape[1],'seconds',round(time.time()-start,1),flush=True)
G=Z.T@I*op.dx;G=(G+G.T)/2;H2=I.T@I*op.dx;H2=(H2+H2.T)/2;del I,Z
src=(root/'ep_kernel_receiving.py').read_text()
np.savez_compressed(root/'EP-kernel-full-reopt-gram.npz',G=G,H2=H2,res=res)
exec(compile(src[src.index("b0=load('DZ','trial31/bs_trial.npz')"):src.index('E0=(1-eps)')], 'penalties', 'exec'))
E0=(1-eps)*q.T@q+x.T@K@q+q.T@K.T@x+x.T@D@x;Ec=cross@W;Ex=K@q+D@x
E=np.zeros_like(G); E[:448,:448]=np.block([[E0,Ec.T,Ex.T],[Ec,selfE,cross[:,119:]],[Ex,cross[:,119:].T,D]])
records=[]
for ii,aa in enumerate((a,ap)):
 delta=float(aa['delta']);t=float(aa['t']);mass=float(aa['new_block_weight']);ws=[(1-mass)*float(w) for w in weights0[:-1]]+[mass,(1-mass)*float(weights0[-1])];rho=np.empty(368)
 for j,b in enumerate(en['physical_blocks']):rho[b['lo']:b['hi']]=(1+1/delta)*float(F(b['allowance']))**2/ws[j+1]/mu
 head=(1+1/delta)*(float(aq)**2/ws[0]+float(er)**2/ws[-1])/mu;P=E-(G+delta*(1+t)*H2)/mu;P[:64,:64]-=x.T@(rho[:,None]*x)+head*q.T@q;P[80:448,:64]-=rho[:,None]*x;P[:64,80:448]-=(rho[:,None]*x).T;P[80:448,80:448]-=np.diag(rho);P=(P+P.T)/2
 jj=np.arange(ii*32,(ii+1)*32);old=np.load(root/f'EP-kernel-joint-reopt-block-{ii}.npz');ix=np.r_[np.arange(64,80),80+old['existing_indices']];A0=-P[np.ix_(ix,ix)];B0=P[np.ix_(ix,jj)];C0=P[np.ix_(jj,jj)];checks={k:float(np.max(np.abs(v-old[k]))) for k,v in [('A',A0),('B',B0),('C',C0)]};assert checks['A']<1e-7 and checks['B']<1e-9 and checks['C']<1e-9,checks
 rec={'variant':ii,'delta':delta,'t':t,'previous_block_replay_errors':checks,'stages':[]}
 for end in (448,len(P)):
  A=-P[64:end,64:end];B=P[64:end,jj];ev=np.linalg.eigvalsh(A);assert ev[0]>0,ev[0]
  for reserve in (0.,.1,1.,10.):
   Ap=A.copy();Ap[:16,:16]+=reserve*np.eye(16);T=np.linalg.solve(Ap,B);S0=C0+B.T@T;S0=(S0+S0.T)/2;vv=np.linalg.eigvalsh(S0);strong=S0[:-1,:-1];c=S0[:-1,-1];r={'optimized_coordinates':end-64,'reserve':reserve,'hessian_min':float(ev[0]),'min_eig':float(vv[0]),'second_min_eig':float(vv[1]),'positive_count':int(sum(vv>0)),'frozen_ray':float(S0[-1,-1]),'schur':float(S0[-1,-1]-c@np.linalg.solve(strong,c)),'kernel_coef_norm':float(np.linalg.norm(T[:16],2))};rec['stages'].append(r);print(ii,r,flush=True)
   if reserve==0:np.savez_compressed(root/f'EP-kernel-full-reopt-block-{ii}-{end-64}.npz',A=A,B=B,C=C0,S=S0,T=T,Q=aa['Q'],eigenvalues=vv)
 records.append(rec)
out={'scope':'Frozen 32-dimensional EO head trials, joint numerical optimization of all 368 inherited physical coordinates, 16 additional kernel coordinates, and optionally 116 reconstructed exceptional coordinates; inherited penalties retained, new kernel action/energy/tail enclosures and their penalty mass omitted. Artificial diagonal reserves are diagnostics, not certified bounds.','sample_count':op.N,'exceptional_rank':rank,'solves':len(res),'max_relative_solve_residual':max(res),'records':records};(root/'EP-kernel-full-reopt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
