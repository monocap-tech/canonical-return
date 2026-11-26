from pathlib import Path
import numpy as np
from fractions import Fraction
from scipy.sparse.linalg import LinearOperator,cg
from scipy.linalg import cho_factor,cho_solve
P=Path('..'); S=10**12;N=32768;d=119;m=192
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
print('new eig',np.linalg.eigvalsh(M)[-24:],'old Xnorm',np.linalg.norm(oldX@q),'new',np.linalg.norm(X),flush=True)
# represent X as full 16x119 trial for exact multiplication by rational Q, arbitrary projection coefficients in raw F and g basis
ev,Q=np.linalg.eigh(M);ntrial=19;q=np.rint(Q[:,-ntrial:][:,::-1]*S).astype(np.int64);qf=q/S; print('selected',ev[-ntrial:],flush=True);V=V@qf;Xfull=X;Xgrid=np.rint(Xfull*S).astype(np.int64);vgrid=np.rint(V[:-1]*S).astype(np.int64);ggrid=np.rint(V[-1:]*S).astype(np.int64)
# numerical joint inverse of actual rounded proposed physical representative, derive t and u
B0=np.vstack([qf,Xgrid/S@qf]);a=np.vstack([B0,-vgrid/S,-ggrid/S]);UU=U@a;full=np.zeros((N,ntrial));full[mask]=UU;uc=ana(full);uu=(uc[pos]*np.sqrt(2)*1j).real;t=np.rint(np.sqrt(alpha)[:,None]/(1-alpha)[:,None]*uu[rp]*S).astype(np.int64);idx=np.sort(np.argsort(np.sum(uu*uu,axis=1))[-8192:]);np.savez_compressed('bw_full_proposal.npz',Q=q,X=Xgrid,v0=vgrid,g0=ggrid,t=t,alpha=C['weights'],j=2*idx+1,T=np.rint(uu[idx]*S).astype(np.int64),G=G,central=M,scale=S)
