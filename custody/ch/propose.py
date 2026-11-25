"""Heuristic only: propose a wider interior inverse on the frozen physical v."""
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
import sys
root=Path(sys.argv[1]); S=10**12; N=65536
v=np.load(root/'bu_trial.npz'); b=np.load(root.parent/'cr_int_t1d_bu'/'bu_symbol.npz'); old=np.load(root/'even_joint_proposal.npz')
x=-4+(np.arange(N)+.5)*8/N; mask=np.abs(x)<2;k=(np.fft.fftfreq(N)*N).astype(int); pos=np.arange(N//2);shift=np.exp(1j*np.pi*k*(-4+4/N)/4)
def synth(c): return (np.fft.ifft(c*shift[:,None],axis=0)*N/np.sqrt(8)).real
def ana(y): return np.fft.fft(y,axis=0)*np.sqrt(8)/N/shift[:,None]
def cosine(j,C):
 out=np.zeros((N,C.shape[1]),complex)
 for kk in (j//2,-j//2):
  idx=kk%N; fac=(-1.)**((kk-np.where(idx>=N//2,idx-N,idx))//N);np.add.at(out,idx,C*(fac[:,None]/np.sqrt(2)))
 out[0]+=(1-np.sqrt(2))*C[j==0].sum(0)
 return synth(out)[mask]
V=cosine(v['js'],v['c']/S)+np.cosh(x[mask]/2)[:,None]*(v['gc']/S)
# Mid-symbol continuation is proposal-only, with clamping beyond the known prefix.
bs=1.2*(7.98-(b['au'].astype(float)+b['al'].astype(float))/(2*float(b['scale'])));bs=bs[np.minimum(np.abs(k),len(bs)-1)]
def mv(z):
 f=np.zeros((N,1));f[mask,0]=z
 return synth(ana(f)*(1+bs[:,None]))[mask,0]
L=LinearOperator((mask.sum(),)*2,matvec=mv,dtype=float)
initial=cosine(old['j'],old['T']/S);U=np.empty_like(V)
for i in range(V.shape[1]):
 y,info=gmres(L,V[:,i],x0=initial[:,i],rtol=1e-9,atol=0,restart=50,maxiter=20);res=np.linalg.norm(mv(y)-V[:,i])/np.linalg.norm(V[:,i]);print('proposal',i,'info',info,'relative residual',res,flush=True);assert info==0 and np.isfinite(res);U[:,i]=y
f=np.zeros((N,V.shape[1]));f[mask]=U;c=ana(f)[pos];u=np.sqrt(2)*c.real;u[0]=c[0].real
# Select from the existing 16384 symbol rows; proposal never supplies a proof bound.
u=u[:len(b['ids'])];assert np.array_equal(b['ids'],2*np.arange(len(u)))
idx=np.sort(np.argsort(np.sum(u*u,axis=1))[-8192:]);T=np.rint(u[idx]*S);assert np.isfinite(T).all() and np.max(np.abs(T))<2**62
np.savez_compressed(root/'even_joint_proposal.npz',j=2*idx,T=T.astype(np.int64),scale=S)
print('heuristic 8192 interior proposal saved',flush=True)
