from pathlib import Path
import numpy as np,zipfile,io,hashlib,json
from scipy.special import digamma,polygamma
from scipy.linalg import null_space,qr
from fractions import Fraction as F
from math import isqrt
import mpmath as mp
root=Path(__file__).resolve().parent
pins={'EN':'d84732c72265bfa30f86463f8108c9e39f62880e','EG':'06805c7f3bb20286f42ef00627937990cb9bf308'}
def load(s,n):
 with zipfile.ZipFile(root/f'EP-parent-{s}-evidence.zip') as z:
  m=z.read('sha256.txt');assert hashlib.sha1(b'blob '+str(len(m)).encode()+b'\0'+m).hexdigest()==pins[s]
  h={p.removeprefix('./'):d for d,p in (line.split(maxsplit=1) for line in m.decode().splitlines())};raw=z.read(n);assert hashlib.sha256(raw).hexdigest()==h[n]
  return dict(np.load(io.BytesIO(raw),allow_pickle=False))
a=load('EN','new/odd_corrector.npz');b=load('EN','merged/odd_merged_corrector.npz');head=load('EG','odd_head_full.npz');check=load('EN','new/odd_U_full.npz')
S=int(a['coefficient_scale']);P=a['P'].astype(float)/S;ids=a['indices'][:139];new=a['new'];mask=np.isin(b['indices'],ids);U=b['coefficient'].astype(float)/float(b['coefficient_scale']);outside=U[~mask];assert np.array_equal(outside.T@outside,np.eye(368));assert np.all(np.count_nonzero(outside,axis=0)==1)
# Orthogonal nullspace proposal, then integer projection using certified right inverse.
Q=null_space(P);assert Q.shape==(139,16);obj=lambda a:np.vectorize(int,otypes=[object])(a);P0=obj(a['P']);Z=obj(a['Z']);zs=int(a['Zscale']);T=np.vectorize(lambda v:int(F(float(v))*S),otypes=[object])(Q)
for _ in range(2):T-=Z@(P0@T)//(S*zs)
def norm(v,s):
 q=sum(int(x)**2 for x in v.flat);r=isqrt(q);return F(r+(r*r<q),s)
pe=F(str(a['Perror']));rho=F(str(a['rho']));Rbound=norm(Z,zs)/(1-rho);e=Rbound*(norm(P0@T,S*S)+pe*norm(T,S));Q=T.astype(float)/S
# Reconstruct numerical raw signed action exactly as archived bq_action source, using digamma and high precision trigamma on selected diagonal modes.
pp=[(n,p) for p in range(2,55) if all(p%i for i in range(2,isqrt(p)+1)) for n in [p**k for k in range(1,7) if p**k<=54]];assert len(pp)==24
wt=np.array([np.log(p)/np.sqrt(n) for n,p in pp]);t=np.log([n for n,p in pp]);mp.mp.dps=35
ns=np.r_[1,ids,new];w=ns*np.pi/4
z=.25+1j*w/2;ps=digamma(z);tri=np.array([float(mp.re(mp.polygamma(1,mp.mpc(.25,float(v)/2)))) for v in w]);
# psi(n)[3] convention is real trigamma, verify against inherited action below.
d=w*ps.imag/2+w*(np.sin(w[:,None]*t)@wt);bd=ps.imag/(2*w)+tri/4;ad=((4-t)*np.cos(w[:,None]*t)-np.sin(w[:,None]*t)/w[:,None])@wt;diag=(bd-ad+2*(ps.real-np.log(np.pi)))/np.sqrt(2)
ex=np.arange(1,40)*2+.5;coef=2*ex**2*np.exp(-2*(4*np.arange(1,40)+1));rr=(.25+w[0]**2)/(.25+w[1:]**2);tr=1/(ex[:,None]**2+w[None,1:]**2)-rr[None,:]/(ex[:,None]**2+w[0]**2)
def action(ms):
 wm=ms*np.pi/4;dm=wm*digamma(.25+1j*wm/2).imag/2+wm*(np.sin(wm[:,None]*t)@wt);den=wm[:,None]**2-w[None,:]**2;equal=ms[:,None]==ns[None,:];raw=2*(dm[:,None]-d[None,:])/np.where(equal,1,den)/np.sqrt(2);raw=np.where(equal,diag[None,:],raw);return raw[:,1:]-raw[:,:1]*rr+(coef[None,:]/(ex[None,:]**2+wm[:,None]**2)/np.sqrt(2))@tr
ms=head['indices'];C=a['coefficient'].astype(float)/S;A=action(ms);B=A@C;Ycheck=check['Y'].astype(float)/float(check['scale']);assert np.array_equal(check['indices'],ms[:len(Ycheck)]);err=np.max(abs(B[:len(Ycheck)]-Ycheck));print('control error',err,flush=True)
# Compare original receiver kernel actions to full head span. Exact rank statement for physical correctors only, numerical head containment.
Fh=head['Y'].astype(float)/float(head['scale']);V=A[:,:139]@Q;coefh=np.linalg.lstsq(Fh,V,rcond=None)[0];res=V-Fh@coefh
print('head residual',np.linalg.norm(res),np.linalg.norm(V),np.linalg.svd(res,compute_uv=False),flush=True)
out={'scope':'Exact receiver-kernel correction bound conditional on archived enclosures; numerical signed-action head containment diagnostic only. No new positivity certificate.','receiver_nullity':16,'outside_receiver_permutation_exact':True,'kernel_center_residual':str(norm(P0@T,S*S)),'true_kernel_coefficient_error_upper':str(e),'true_kernel_coefficient_error_upper_float':float(e),'center_kernel_singular_min':float(np.linalg.svd(Q,compute_uv=False)[-1]),'control_max_abs_action_difference':float(err),'control_relative_frobenius':float(np.linalg.norm(B[:len(Ycheck)]-Ycheck)/np.linalg.norm(Ycheck)),'head_rank':int(np.linalg.matrix_rank(Fh)),'kernel_action_frobenius':float(np.linalg.norm(V)),'head_fit_residual_frobenius':float(np.linalg.norm(res)),'relative_head_fit_residual':float(np.linalg.norm(res)/np.linalg.norm(V)),'head_fit_residual_singular_values':np.linalg.svd(res,compute_uv=False).tolist(),'kernel_action_singular_values':np.linalg.svd(V,compute_uv=False).tolist(),'head_coeff_norm':float(np.linalg.norm(coefh))}
(root/'EP-receiver-kernel-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
# Test against the complete existing action bank, including all 368 physical columns.
with zipfile.ZipFile(root/'EP-parent-EG-evidence.zip') as z: eg=json.loads(z.read('extension.json'))
with zipfile.ZipFile(root/'EP-parent-EN-evidence.zip') as z: en=json.loads(z.read('extension.json'))
prefix=load('EN','merged/odd_combined_projection.npz')['Y']; old=[prefix[:,119:471]];new=[]
for c in eg['new_chunks']:
 ar=load('EG',f"odd_U_{c['lo']}_{c['hi']}.npz");old.append(ar['Y'])
for c in en['new_chunks']:
 ar=load('EN',f"new/odd_U_{c['lo']}_{c['hi']}.npz");new.append(ar['Y'])
Y=np.column_stack([Fh,np.vstack(old)/1e12,np.vstack(new)/1e12]);assert Y.shape==(32768,487);assert np.array_equal(np.rint(Y[:8192]*1e12).astype(np.int64),prefix)
# Thin QR avoids an enormous full left SVD.
qb,rb=np.linalg.qr(Y,mode='reduced');sv=np.linalg.svd(rb,compute_uv=False);res2=V-qb@(qb.T@V);s2=np.linalg.svd(res2,compute_uv=False)
gramdef=norm(T.T@T-np.eye(16,dtype=object)*S*S,S*S);assert gramdef< F(1,100);assert e<F(1,100)
# True basis W=(I-R P*) T/S; sigma_min(W)>=sqrt(1-gramdef)-e>0.98.
assert gramdef<F(1,100) and e<F(1,100)
out.update(combined_action_rank=int(np.sum(sv>sv[0]*max(Y.shape)*np.finfo(float).eps)),combined_action_singular_min=float(sv[-1]),combined_fit_residual_frobenius=float(np.linalg.norm(res2)),combined_fit_relative=float(np.linalg.norm(res2)/np.linalg.norm(V)),combined_fit_residual_singular_values=s2.tolist(),center_gram_defect_exact=str(gramdef),center_gram_defect_float=float(gramdef),true_kernel_basis_sigma_min_lower='49/50',true_kernel_basis_rank=16)
(root/'EP-receiver-kernel-audit.json').write_text(json.dumps(out,indent=2)+'\n')
np.savez_compressed(root/'EP-receiver-kernel-proposal.npz',indices=ids,coefficient=np.array(T,str),coefficient_scale=str(S),coefficient_error=str(e),action_center=V,action_indices=ms)
print('combined', {k:v for k,v in out.items() if k.startswith('combined') or k.endswith('_float') or k.startswith('true_kernel_basis')},flush=True)
