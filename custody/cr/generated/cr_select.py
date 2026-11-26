from pathlib import Path
import numpy as np
from fractions import Fraction
from scipy.sparse.linalg import LinearOperator,cg
from scipy.linalg import cho_factor,cho_solve
P=Path('..'); S=10**12;N=32768;d=119;m=112
A=np.load(Path('odd_combined_projection.npz'));Y=A['Y']/S;B=np.load(P/'cr_int_t1d_bs/bs_trial.npz');q=np.eye(d);oldX=B['X']/S;E=np.load(Path('odd_BQ_energy.npz'));K=E['K'].astype(float)/float(E['scale']);D=E['D'].astype(float)/float(E['scale']);mu=float(Fraction(str(B['mu'])));eps=float(Fraction(str(B['eps'])))
x=-4+(np.arange(N)+.5)*8/N;mask=abs(x)<2;xx=x[mask];k=np.fft.fftfreq(N)*N;k=k.astype(int);pos=np.arange(N//2);neg=(-pos-1)%N;phase=np.exp(1j*np.pi*.5*x/4);shift=np.exp(1j*np.pi*k*(-4+4/N)/4)
def synth(C):return (np.fft.ifft(C*shift[:,None],axis=0)*N*phase[:,None]/np.sqrt(8)).real
def ana(Y):return np.fft.fft(Y/phase[:,None],axis=0)*np.sqrt(8)/N/shift[:,None]
# arbitrary circular sine synth, including aliased high dictionary frequencies (proposal only)
def sin_synth(j,C):
 out=np.zeros((N,C.shape[1]),complex);k0=(j-1)//2;k1=(-j-1)//2
 for kk,sgn in [(k0,1),(k1,-1)]:
  idx=kk%N;fac=(-1.)**((kk- np.where(idx>=N//2,idx-N,idx))//N);np.add.at(out,idx,C*(sgn*fac[:,None]/(np.sqrt(2)*1j)))
 return synth(out)[mask]
ms=B['ms'];nc=np.column_stack([Y[:,:d]@q,Y[:,d:]])* (np.sqrt(2)*np.where(((ms+1)//2)%2,-1,1))[:,None];cc=np.zeros((N,nc.shape[1]),complex);cc[ms]=nc/(np.sqrt(2)*1j);cc[(-ms)%N]=-nc/(np.sqrt(2)*1j);Ys=(np.fft.ifft(cc*shift[:,None],axis=0)*N/np.sqrt(8)).real[mask]
PD=np.load(P/'cr_int_t1d_bk/odd/physical_dictionary.npz');EE=np.load(P/'cr_int_t1d_bk/odd/exceptional.npz');Fraw=sin_synth(PD['indices'], np.sqrt(PD['weights']/float(PD['weight_scale']))[:,None]*(EE['E']/float(EE['scale'])));g=np.sinh(xx/2)[:,None];Z=np.column_stack([Ys,Fraw,g]);print('sample shape',Z.shape,flush=True); assert np.max(abs(Ys.T@Ys*8/N - np.column_stack([Y[:,:d]@q,Y[:,d:]]).T@np.column_stack([Y[:,:d]@q,Y[:,d:]])))<1e-8
C=np.load(P/'cr_int_t1d_br/odd_surplus_cross.npz');rj=C['indices'];alpha=C['weights']/1e8;rp=rj//2;rn=(-rp-1)%N
b=np.load(P/'cr_int_t1d_bs/bs_symbol.npz');bvals=1.2*(7.98-(b['au'].astype(float)+b['al'].astype(float))/(2*float(b['scale'])));bs=bvals[np.minimum(np.abs(2*k+1)//2,len(bvals)-1)];bs[rp]+=alpha/(1-alpha);bs[rn]+=alpha/(1-alpha)
def mv(v):
 a=np.zeros((N,1));a[mask,0]=v;return synth(ana(a)*(1+bs[:,None]))[mask,0]
L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float);U=np.empty_like(Z)
for i in range(Z.shape[1]):
 v,info=cg(L,Z[:,i],rtol=1e-9,maxiter=160);assert info==0;U[:,i]=v
 if i%24==0:print('solved',i,flush=True)
G=Z.T@U*8/N;G=(G+G.T)/2;n=d+m;proj=np.linalg.solve(G[n:,n:],G[n:,:n]);J=G[:n,:n]-G[:n,n:]@proj;J=(J+J.T)/2
# maximize lower bracket energy - weighted residual /mu; small inflation for numerical safety
fac=1.0;X=np.linalg.solve(fac*J[d:,d:]/mu-D,K@q-fac*J[d:,:d]/mu);V=proj@np.vstack([np.eye(d),X]);E0=(1-eps)*q.T@q+q.T@K.T@X+X.T@K@q+X.T@D@X;M=E0-fac*np.vstack([np.eye(d),X]).T@J@np.vstack([np.eye(d),X])/mu
print('new eig',np.linalg.eigvalsh(M)[-16:],'old Xnorm',np.linalg.norm(oldX@q),'new',np.linalg.norm(X),flush=True)

ev,Q=np.linalg.eigh(M);q=Q[:,-14:-13];a=np.vstack([q,X@q,-V@q]);u=U@a
old=np.load('odd_merged_corrector.npz')['indices']
cand=np.arange(3,4096,2);cand=cand[~np.isin(cand,old)];score=np.empty(len(cand))
for i in range(0,len(cand),128):
    ns=cand[i:i+128];basis=np.sin(np.pi*xx[:,None]*ns[None,:]/4)/np.sqrt(2)
    score[i:i+len(ns)]=np.linalg.norm(basis.T@u*8/N,axis=1)*np.log(2+ns)**.5
order=np.argsort(score)[-16:][::-1];new=cand[order]
assert len(old)==251 and len(new)==16 and not np.intersect1d(new,old).size
print('CR targets',ev[-14:-13].tolist(), 'selected',new.tolist(), 'scores',score[order].tolist(),flush=True)
np.save('odd_new.npy',new);np.savez_compressed('cr_selection.npz',cand=cand,score=score,eig=ev[-14:-13],new=new,direction=q)
