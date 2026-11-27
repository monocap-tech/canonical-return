from pathlib import Path
import zipfile,hashlib,io,json
from fractions import Fraction as F
from math import isqrt
import numpy as np
root=Path(__file__).resolve().parent
pins={'EO':'7e519ccdca9440b65103fffa20b24d8e6edeba9b','EN':'d84732c72265bfa30f86463f8108c9e39f62880e','EG':'06805c7f3bb20286f42ef00627937990cb9bf308','DV':'dc58e2c0003c4d47ccaec7ab2251a6af4890f48b','DZ':'036d47d3710631055e0cc6e7458af7ab339ae688','CR':'85b77e4beaddbeba6c69a40087db283bc52c182f','DA':'5bd5ba86d8454e6c730c03f4ed890ef3587957f1','DM':'843341d57963dd0f37da1f6eaefbc5e6bd537dad'}
archives={}; z={}
for stage,pin in pins.items():
 p=root/f'EP-parent-{stage}-evidence.zip'; z[stage]=zipfile.ZipFile(p); a=z[stage]; assert a.testzip() is None
 m=a.read('sha256.txt'); assert hashlib.sha1(b'blob '+str(len(m)).encode()+b'\0'+m).hexdigest()==pin
 rows=m.decode().splitlines()
 for row in rows:
  h,n=row.split(maxsplit=1);assert hashlib.sha256(a.read(n.removeprefix('./'))).hexdigest()==h
 archives[stage]={'bytes':p.stat().st_size,'entries':len(rows),'manifest_blob':pin}
def load(stage,path):
 with np.load(io.BytesIO(z[stage].read(path)),allow_pickle=False) as a:return {k:a[k] for k in a.files}
def report(stage,path):return json.loads(z[stage].read(path))
d=119;m=352;add=16;ids=np.arange(1,65536,2); eg=report('EG','extension.json');en=report('EN','extension.json')
h=load('EG','odd_head_full.npz');dv=load('DV','merged/odd_combined_projection.npz');e=load('EN','merged/odd_combined_projection.npz')
assert np.array_equal(h['indices'],ids) and h['Y'].shape==(32768,d)
assert np.array_equal(h['Y'][:8192],dv['Y'][:,:d]);assert int(h['scale'])==int(dv['Yscale'])==int(e['Yscale'])==10**12
old_chunks=[];new_chunks=[];old_ids=[];new_ids=[]
for b in eg['new_chunks']:
 lo,hi=b['lo'],b['hi'];a=load('EG',f'odd_U_{lo}_{hi}.npz');assert np.array_equal(a['indices'],np.arange(lo,hi+1,2));assert a['Y'].shape==(len(a['indices']),m) and int(a['scale'])==10**12;old_chunks.append(a['Y']);old_ids.extend(a['indices'])
for b in en['new_chunks']:
 lo,hi=b['lo'],b['hi'];a=load('EN',f'new/odd_U_{lo}_{hi}.npz');assert np.array_equal(a['indices'],np.arange(lo,hi+1,2));assert a['Y'].shape==(len(a['indices']),add) and int(a['scale'])==10**12;new_chunks.append(a['Y']);new_ids.extend(a['indices'])
assert np.array_equal(old_ids,ids[8192:]) and np.array_equal(new_ids,ids)
old=np.vstack([dv['Y'][:,d:],*old_chunks]);new=np.vstack(new_chunks);del old_chunks,new_chunks
assert np.array_equal(np.column_stack([h['Y'][:8192],old[:8192],new[:8192]]),e['Y'])
assert en['physical_blocks'][:-1]==eg['physical_blocks']
old_e=load('DV','merged/odd_BQ_energy.npz');new_e=load('EN','merged/odd_BQ_energy.npz')
# Compare exact rational centers in Frobenius norm against sum of stored global allowances.
energy={}
for key in ['K','D']:
 a=old_e[key];b=new_e[key][:m,:a.shape[1]];s0=int(old_e['scale']);s1=int(new_e['scale']);scale=s0*s1
 sq=sum((int(x)*s1-int(y)*s0)**2 for x,y in zip(a.flat,b.flat));r=isqrt(sq);upper=F(r+(r*r<sq),scale);allow=F(str(old_e[key+'error']))+F(str(new_e[key+'error']));assert upper<=allow
 energy[key]={'center_difference_norm_upper':float(upper),'sum_global_error_allowances':float(allow),'exact_upper_le_allowance':True}
from scipy.sparse.linalg import LinearOperator,cg
S=10**12; dv,dz=[load(s,'trial31/bs_trial.npz') for s in ('DV','DZ')]
assert np.array_equal(dv['js'],dz['js']) and np.array_equal(dv['a'],dz['a'])
rj=dv['js'][-len(dv['a']):];alpha=dv['a']/1e8;jf=dv['js'][:-len(rj)]
assert len(rj)==1882 and len(jf)==5143 and len(np.unique(dv['js']))==7025 and np.all(rj%2==1) and np.all((alpha>0)&(alpha<1))
aroot=np.array([isqrt(int(v)*S*S//10**8) for v in dv['a']],object)
for a in (dv,dz):assert np.array_equal(a['c'][-len(rj):],-(aroot[:,None]*a['t'].astype(object)//S).astype(np.int64))
class Op:
 def __init__(self,N):
  self.N=N;self.x=-4+(np.arange(N)+.5)*8/N;self.mask=abs(self.x)<2;self.xx=self.x[self.mask];self.dx=8/N
  k=(np.fft.fftfreq(N)*N).astype(int);self.shift=np.exp(1j*np.pi*k*(-4+4/N)/4);self.phase=np.exp(1j*np.pi*.5*self.x/4)
  j=2*np.arange(N//2)+1;pp=[]
  for p in range(2,55):
   if all(p%i for i in range(2,int(np.sqrt(p))+1)):
    a=p
    while a<=54:pp.append((a,p));a*=p
  assert len(pp)==24
  av=sum(2*np.log(p)/np.sqrt(a)*np.cos(np.pi*j*np.log(a)/8) for a,p in pp)
  bs=1.2*(7.98-av[np.abs(2*k+1)//2]);rp=rj//2;rn=(-rp-1)%N;bs[rp]+=alpha/(1-alpha);bs[rn]+=alpha/(1-alpha);self.mult=1+bs
  print("full multiplier min",float(self.mult.min()),flush=True)
  self.L=LinearOperator((len(self.xx),)*2,matvec=self.mv,dtype=float)
 def mv(self,v):
  a=np.zeros(self.N);a[self.mask]=v
  return (np.fft.ifft(np.fft.fft(a/self.phase)*self.mult)*self.phase).real[self.mask]
 def sin(self,j,C):
  out=np.empty((len(self.xx),C.shape[1]));N=self.N
  for b in range(0,C.shape[1],8):
   cc=C[:,b:b+8];v=np.zeros((N,cc.shape[1]),complex)
   for kk,sgn in [((j-1)//2,1),((-j-1)//2,-1)]:
    idx=kk%N;fac=(-1.)**((kk-np.where(idx>=N//2,idx-N,idx))//N)
    np.add.at(v,idx,cc*(sgn*fac[:,None]/(np.sqrt(2)*1j)))
   out[:,b:b+8]=(np.fft.ifft(v*self.shift[:,None],axis=0)*N*self.phase[:,None]/np.sqrt(8)).real[self.mask]
  return out
 def neu(self,ms,C):
  N=self.N;v=np.zeros((N,C.shape[1]),complex);nc=C*(np.sqrt(2)*np.where(((ms+1)//2)%2,-1,1))[:,None]
  v[ms]=nc/(np.sqrt(2)*1j);v[(-ms)%N]=-nc/(np.sqrt(2)*1j)
  return (np.fft.ifft(v*self.shift[:,None],axis=0)*N/np.sqrt(8)).real[self.mask]
 def solve(self,b):
  x,info=cg(self.L,b,rtol=2e-11,maxiter=300);r=np.linalg.norm(self.mv(x)-b)/np.linalg.norm(b);assert info==0 and r<1e-9,(info,r);return x,float(r)
op=Op(32768);p=load('DZ','trial31/bw_full_proposal.npz');Y=load('DV','merged/odd_combined_projection.npz')['Y']/S
ix=[0,15,30];q=dz['Q'][:,ix]/S;x=dz['X']/S@q;v=dz['v0'][:,ix]/S;g=dz['g0'][:,ix]/S
W=np.vstack([q,x,-v,-g]);expected=W.T@p['G']@W
B=op.neu(dz['ms'],Y@np.vstack([q,x]))+op.sin(jf,dz['c'][:-len(rj),ix]/S)-np.sinh(op.xx/2)[:,None]*g
U,rs=zip(*(op.solve(b) for b in B.T));actual=B.T@np.column_stack(U)*op.dx;err=float(np.max(abs(actual-expected)));rel=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected));print('historical',err,rel,flush=True);assert rel<1e-6,(err,rel)
trials=[]
for stage,path in [('DV','trial31'),('DV','trial32'),('DZ','trial31'),('CR','trial14'),('DA','trial19')]:
 a=load(stage,path+'/bs_trial.npz')
 assert np.array_equal(a['js'],dv['js']) and np.array_equal(a['a'],dv['a']) and int(a['scale'])==S
 assert np.array_equal(a['c'][-len(rj):],-(aroot[:,None]*a['t'].astype(object)//S).astype(np.int64))
 trials.append(a)
op2=Op(131072);features=op2.sin(jf,-np.column_stack([a['c'][:-len(rj)] for a in trials])/S);features=np.column_stack([features,np.sinh(op2.xx/2)]);u,s,vh=np.linalg.svd(features,full_matrices=False);keep=s>s[0]*1e-7;E=u[:,keep];print('rank',sum(keep),flush=True)
sol=[];rr=[]
for i,b in enumerate(E.T):
 a,r=op2.solve(b);sol.append(a);rr.append(r)
 if i%10==0:print('solves',i,flush=True)
IE=np.column_stack(sol);T=E.T@IE*op2.dx;T=(T+T.T)/2;assert np.linalg.eigvalsh(T)[0]>0
Ya=np.column_stack([h['Y'],old,new])/S;K=new_e['K'].astype(float)/float(new_e['scale']);D=new_e['D'].astype(float)/float(new_e['scale']);mu=float(F(str(dz['mu'])));eps=float(F(str(dz['eps'])))
result={'scope':'Numerical sampled multiplier replay; expanded nested rounded exceptional span only; no directed signed-symbol interval replay or positivity claim','archives':archives,'surplus_indices':len(rj),'dictionary_indices':len(jf),'historical':{'N':32768,'columns':ix,'gram_max_abs_error':err,'gram_relative_error':rel,'max_solve_residual':max(rs),'expected':expected.tolist(),'actual':actual.tolist(),'multiplier_min':float(op.mult.min()),'multiplier_max':float(op.mult.max())},'extended':{'N':131072,'exceptional_columns':features.shape[1],'retained_rank':int(sum(keep)),'relative_svd_cutoff':1e-7,'max_solve_residual':max(rr),'multiplier_min':float(op2.mult.min()),'multiplier_max':float(op2.mult.max()),'variants':{}}}
for name in ('candidate32.npz','candidate32_positive_t.npz'):
 a=load('EO',name);q=a['Q'][:,-1];x=a['X']@q;bb=op2.neu(h['indices'],(Ya@np.r_[q,x])[:,None])[:,0];ib,r=op2.solve(bb);gamma=float(a['delta'])*(1+float(a['t']));bgram=float(bb@ib*op2.dx);hh=E.T@ib*op2.dx
 def evaluate(gamma,rank):
  Tg=T[:rank,:rank]+gamma*(IE[:,:rank].T@IE[:,:rank])*op2.dx;hg=hh[:rank]+gamma*(IE[:,:rank].T@ib)*op2.dx;bg=bgram+gamma*(ib@ib)*op2.dx;c=np.linalg.solve(Tg,hg);metric=bg-hg@c;physical=(1-eps)*(q@q)+2*x@K@q+x@D@x
  return {'unprojected_metric':float(bg),'restricted_projected_metric':float(metric),'physical_energy':float(physical),'energy_minus_metric_over_mu':float(physical-metric/mu),'projection_gain':float(hg@c),'gamma':gamma}
 result['extended']['variants'][name]={'solve_relative_residual':r,'nested':[{'cutoff':cut,'rank':int(sum(s>s[0]*cut)), 'ideal':evaluate(0,int(sum(s>s[0]*cut))), 'gamma_metric':evaluate(gamma,int(sum(s>s[0]*cut)))} for cut in (1e-5,1e-6,1e-7)]}
(root/'EP-expanded-exceptional-replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
