"""EK: higher-resolution numerical proposal screen, not a certificate.
Run only after EJ acceptance, with completed EI evidence attached and CE recovered.
All sampled FFT, CG and eigenvalue calculations are proposal-only.
"""
from pathlib import Path
import os,json
from fractions import Fraction as F
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg
root=Path(__file__).resolve().parents[2];w=root/'ce-output/pristine';cb=w/'cr_int_t1d_cb';eg=root/'custody/eg/generated';out=root/'custody/ek/generated';out.mkdir(parents=True,exist_ok=True)
ext=json.loads((eg/'extension.json').read_text())
audit=json.loads((root/'custody/ej/generated/audit.json').read_text());assert audit['best_margin_audited']
ei=json.loads((root/'custody/ei/generated/result.json').read_text());assert ei['best']['positive'] and ei['dimension']==31
assert ext['cutoff']==65536 and ext['head_half_target_met'] and ext['all_block_half_targets_met']
S=10**12;N=131072;d=119;m=352;n=d+m
old=np.load(root/'custody/dv/generated/merged/odd_combined_projection.npz')
h=np.load(eg/'odd_head_full.npz');ms=h['indices']
assert np.array_equal(ms,np.arange(1,65536,2))
assert h['Y'].shape==(32768,d) and np.array_equal(h['Y'][:8192],old['Y'][:,:d])
chunks=[];covered=[]
for b in ext['new_chunks']:
    lo,hi=b['lo'],b['hi'];a=np.load(eg/f'odd_U_{lo}_{hi}.npz')
    assert np.array_equal(a['indices'],np.arange(lo,hi+1,2)) and a['Y'].shape==(len(a['indices']),m)
    chunks.append(a['Y']);covered.extend(map(int,a['indices']))
assert covered==list(range(16385,65536,2))
Y=np.column_stack([h['Y'],np.vstack([old['Y'][:,d:],*chunks])]).astype(float)/S
del chunks
E=np.load(cb/'odd_BQ_energy.npz');K=E['K'].astype(float)/float(E['scale']);D=E['D'].astype(float)/float(E['scale'])
B=np.load(cb/'bs_trial.npz');mu=float(F(str(B['mu'])));eps=float(F(str(B['eps'])))
x=-4+(np.arange(N)+.5)*8/N;mask=abs(x)<2;xx=x[mask]
k=(np.fft.fftfreq(N)*N).astype(int);pos=np.arange(N//2)
phase=np.exp(1j*np.pi*.5*x/4);shift=np.exp(1j*np.pi*k*(-4+4/N)/4)
def synth(C):return (np.fft.ifft(C*shift[:,None],axis=0)*N*phase[:,None]/np.sqrt(8)).real
def ana(v):return np.fft.fft(v/phase[:,None],axis=0)*np.sqrt(8)/N/shift[:,None]
# Tile FFTs to avoid a multi-gigabyte temporary array.
def sin_synth(j,C):
    result=np.empty((len(xx),C.shape[1]))
    k0=(j-1)//2;k1=(-j-1)//2
    for b in range(0,C.shape[1],16):
        cc=C[:,b:b+16];v=np.zeros((N,cc.shape[1]),complex)
        for kk,sgn in [(k0,1),(k1,-1)]:
            idx=kk%N;fac=(-1.)**((kk-np.where(idx>=N//2,idx-N,idx))//N)
            np.add.at(v,idx,cc*(sgn*fac[:,None]/(np.sqrt(2)*1j)))
        result[:,b:b+16]=synth(v)[mask]
    return result
Ys=np.empty((len(xx),n))
for b in range(0,n,16):
    nc=Y[:,b:b+16]*(np.sqrt(2)*np.where(((ms+1)//2)%2,-1,1))[:,None]
    v=np.zeros((N,nc.shape[1]),complex);v[ms]=nc/(np.sqrt(2)*1j);v[(-ms)%N]=-nc/(np.sqrt(2)*1j)
    Ys[:,b:b+16]=(np.fft.ifft(v*shift[:,None],axis=0)*N/np.sqrt(8)).real[mask]
gram_err=float(np.max(abs(Ys.T@Ys*8/N-Y.T@Y)));assert gram_err<1e-7
PD=np.load(w/'cr_int_t1d_bk/odd/physical_dictionary.npz');EE=np.load(w/'cr_int_t1d_bk/odd/exceptional.npz')
Fraw=sin_synth(PD['indices'],np.sqrt(PD['weights']/float(PD['weight_scale']))[:,None]*(EE['E']/float(EE['scale'])))
Z=np.column_stack([Ys,Fraw,np.sinh(xx/2)]);del Ys,Fraw,Y
C=np.load(w/'cr_int_t1d_br/odd_surplus_cross.npz');rj=C['indices'];alpha=C['weights']/1e8;rp=rj//2;rn=(-rp-1)%N
assert np.all((rp>=0)&(rp<N//2))
# Independent numerical extension of the signed symbol. Match every original
# stored directed interval before using the higher frequencies for proposals.
j=2*pos+1
pp=[]
for p in range(2,55):
    if all(p%i for i in range(2,int(np.sqrt(p))+1)):
        a=p
        while a<=54:pp.append((a,p));a*=p
assert len(pp)==24
av=sum(2*np.log(p)/np.sqrt(a)*np.cos(np.pi*j*np.log(a)/8) for a,p in pp)
oldb=np.load(w/'cr_int_t1d_bs/bs_symbol.npz');l=oldb['al']/float(oldb['scale']);u=oldb['au']/float(oldb['scale'])
assert np.array_equal(oldb['ids'],j[:len(l)])
symbol_violation=float(max(0,np.max(l-av[:len(l)]),np.max(av[:len(l)]-u)))
assert symbol_violation<1e-8,('signed symbol mismatch',symbol_violation)
bvals=1.2*(7.98-av);bs=bvals[np.abs(2*k+1)//2];bs[rp]+=alpha/(1-alpha);bs[rn]+=alpha/(1-alpha)
def mv(v):
    a=np.zeros((N,1));a[mask,0]=v
    return synth(ana(a)*(1+bs[:,None]))[mask,0]
L=LinearOperator((len(xx),)*2,matvec=mv,dtype=float);U=np.empty_like(Z)
for i in range(Z.shape[1]):
    v,info=cg(L,Z[:,i],rtol=1e-9,maxiter=200);assert info==0,(i,info);U[:,i]=v
    if i%24==0:print('EK solved',i,flush=True)
G=Z.T@U*8/N;G=(G+G.T)/2
Hinv=U.T@U*8/N;Hinv=(Hinv+Hinv.T)/2
del Z
oldp=np.load(root/'custody/dv/generated/proposal32.npz')
assert G.shape==oldp['G'].shape
best=ei['best'];delta0,t0=map(lambda s:float(F(s)),best['parameters'][1:3])
ws=list(map(lambda s:float(F(s)),best['weights']));assert len(ws)==21 and abs(sum(ws)-1)<1e-12
aq=float(F(ei['head_allowance']));er=float(F(ei['representation_error']))
rho0=np.empty(m)
for i,b in enumerate(ext['physical_blocks']):
    rho0[b['lo']:b['hi']]=float(F(b['allowance']))**2/ws[i+1]
def calc(delta,t,ideal=False):
    gamma=0 if ideal else delta*(1+t)
    H=G+gamma*Hinv
    proj=np.linalg.solve(H[n:,n:],H[n:,:n]);J=H[:n,:n]-H[:n,n:]@proj;J=(J+J.T)/2
    rho=np.zeros(m) if ideal else (1+1/delta)*rho0/mu
    head=0 if ideal else (1+1/delta)*(aq*aq/ws[0]+er*er/ws[-1])/mu
    A=J[d:,d:]/mu-D+np.diag(rho);A=(A+A.T)/2
    amin=float(np.linalg.eigvalsh(A)[0]);assert amin>0
    X=np.linalg.solve(A,K-J[d:,:d]/mu);V=proj@np.vstack([np.eye(d),X]);W=np.vstack([np.eye(d),X,-V])
    M=(1-eps)*np.eye(d)+K.T@X+X.T@K+X.T@D@X-W.T@G@W/mu;M=(M+M.T)/2
    P=M-gamma*(W.T@Hinv@W)/mu-X.T@(rho[:,None]*X)-head*np.eye(d);P=(P+P.T)/2
    ev,Q=np.linalg.eigh(P)
    R=Q[:,-32:][:,::-1]
    record=dict(delta=delta,t=t,ideal=ideal,proxy31=float(ev[-31]),proxy32=float(ev[-32]),proxy33=float(ev[-33]),nominal_on_trial32=float(np.linalg.eigvalsh(R.T@M@R)[0]),hessian_min=amin)
    return record,(X,V,R)
ideal,_=calc(1,0,True);rows=[];candidate=None;bestrec=None
for delta in sorted(set([delta0]+[delta0*10**v for v in np.linspace(-1,1,9)])):
    for t in (0,t0,2*t0):
        rec,arr=calc(delta,t);rows.append(rec)
        if bestrec is None or rec['proxy32']>bestrec['proxy32']:bestrec,candidate=rec,arr
X,V,R=candidate
np.savez_compressed(out/'candidate32.npz',X=X,V=V,Q=R,delta=bestrec['delta'],t=bestrec['t'])
report=dict(scope='sampled numerical screen only; omitted exact residual and enclosure terms, t=0 not admissible in exact comparison',physical_cutoff=65536,sample_count=N,gram_replay_error=gram_err,symbol_interval_violation=symbol_violation,old_nominal32=float(np.linalg.eigvalsh(oldp['central'])[-32]),ideal=ideal,best=bestrec,rows=rows)
(out/'screen.json').write_text(json.dumps(report,indent=2)+'\n')
print('EK SPECTRUM',ideal,bestrec,flush=True)
