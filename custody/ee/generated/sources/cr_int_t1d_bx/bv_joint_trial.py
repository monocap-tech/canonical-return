from bv_response import *
# v and u are represented exactly by their saved rational-grid synthesis centers;
# all uncertainty of v vs actual projected physical head is applied to final inverse norm.
def cross(A,Ae,B,Be,w):
 return A.T@(w[:,None]*B),abs(A).T@(abs(w)[:,None]*Be)+Ae.T@(abs(w)[:,None]*abs(B))+Ae.T@(abs(w)[:,None]*Be)
def certificate(d,a,C=None):
 G=d['G'];Ge=rows(d['Ge']);q=d['q'];x=d['x'];diag=d['diag'];Z=d['Z'];Ze=d['Ze'];ZS=d['ZS'];assert ZS==S*S;assert np.array_equal(a['ids'],d['c']['ids']);U=obj(a['U']);Ue=obj(a['Ue']);Gu=cf(obj(a['G']),S**3);Gue=rows(cf(obj(a['Ge']),S**3));T=obj(a['T']);ix=np.searchsorted(a['ids'],a['j']);Gvu=cf(Z[ix].T@T,S**3);Gvue=rows(cf(Ze[ix].T@abs(T),S**3)) # row norm not sufficient for rectangular operator unless use both max row and col or Frobenius
 Gvue=norm(Ze[ix].T@abs(T),S**3)
 Z2u=cf(obj(a['Z2']),S**4);Z2ue=rows(cf(obj(a['Z2e']),S**4));ut=Gu-Z2u;ute=Gue+Z2ue;utn=rows(ut)+ute;assert utn>=0
 # Exact B cross interval using global signed symbol bounds on omitted tail.
 w=np.array([scalarfloor((l+h)/2) for l,h in zip(d['bl'],d['bh'])],object);we=max(max(abs(F(int(v),S)-l),abs(F(int(v),S)-h)) for v,l,h in zip(w,d['bl'],d['bh']));mid=(d['bmin']+d['bmax'])/2;rad=(d['bmax']-d['bmin'])/2
 c,ce=cross(Z,Ze,U,Ue,w);t,te=cross(Z,Ze,U,Ue,np.ones(len(U),object));Vut=Gvu-cf(t,S**4);Vute=Gvue+norm(te,S**4);Bvu=cf(c,S**5)+mid*Vut;Bvue=norm(ce,S**5)+we*sqrtup((rows(d['Z2'])+rows(d['Z2E']))*(rows(Z2u)+Z2ue))+abs(mid)*Vute+rad*sqrtup(d['vt']*utn)
 c,ce=gram(U,Ue,w);Bu=cf(c,S**5)+mid*ut;Bue=rows(cf(ce,S**5))+we*(rows(Z2u)+Z2ue)+abs(mid)*ute+rad*utn
 # physical square upper J² minus exterior lower; all intervals directed.
 w2=np.array([scalarceil(max(l*l,h*h)) for l,h in zip(d['bl'],d['bh'])],object);j2,j2e=gram(U,Ue,w2);b2max=max(d['bmin']**2,d['bmax']**2);J2=cf(j2,S**5)+b2max*ut;J2e=rows(cf(j2e,S**5))+b2max*ute
 HG=cf(obj(a['HG']),S**3);HGe=rows(cf(obj(a['HGE']),S**3));CX=cf(obj(a['CX']),S**3);CE=rows(cf(obj(a['CE'])+obj(a['CE']).T,S**3));Y2=cf(obj(a['Y2']),S**4);Y2e=rows(cf(obj(a['Y2e']),S**4));ht=rows(HG-Y2)+HGe+Y2e;assert ht>=0;O=CX+CX.T-HG;Oe=CE+HGe+2*max(abs(d['bmin']),d['bmax'])*sqrtup(utn*ht);B2u=J2-O;B2ue=J2e+Oe
 print('U tail',float(utn),'Bvu error',float(Bvue),'Bu error',float(Bue),'B2u error',float(B2ue),flush=True)
 if C is None:C=np.rint(np.linalg.solve(np.array(Bu+B2u,float),np.array(Bvu.T,float))*S).astype(np.int64)
 Cq=cf(C,S);nc=norm(C,S);bc=Bvu@Cq;H=G-bc-bc.T+Cq.T@(Bu+B2u)@Cq+diag;eh=Ge+2*nc*Bvue+nc*nc*(Bue+B2ue)
 # Use C=I also for matched direct inverse proposal; compare with weighted optimizing C.
 E=d['E'];dd=d['derr']+d['e'];h=rows(H)+eh;assert h>=0;allow=d['eE']+(eh+2*sqrtup(h)*dd+dd*dd)/d['mu'];L=E-H/d['mu'];low=exact_lower(L)-allow
 print('trial',d['par'],'Cnorm',float(nc),'head',float(h),'allow',float(allow),'lower',float(low),'central eigenvalues',('skipped fixed' if '--fixed' in sys.argv else np.linalg.eigvalsh(np.array(L,float))),flush=True)
 return C,L,allow,low
if __name__=='__main__':
 p=sys.argv[1];d=loadjoint();a=np.load(root/f'{p}_joint_Uaction.npz');fixed='--fixed' in sys.argv
 for tag in (['opt'] if fixed else ['direct','opt']):
  f=root/f'{p}_joint_{tag}_accepted.npz';old=np.load(f) if fixed else None;C=old['C'] if fixed else (np.eye(d['q'].shape[1],dtype=np.int64)*S if tag=='direct' else None);C,L,allow,low=certificate(d,a,C)
  if fixed:R=old['R']
  else:_,R=np.linalg.eigh(np.array(L,float));R=np.rint(R*S).astype(np.int64)
  best=0;bestl=F(0)
  for n in range(1,d['q'].shape[1]+1):
   r=R[:,-n:];rq=cf(r,S);dr=rows(rq.T@rq-I(n));ql=rows(cf(d['q'].T@d['q'],S*S)-I(d['q'].shape[1]));rl=exact_lower(rq.T@L@rq)-(1+dr)*allow;normalized=rl/((1+dr)*(1+ql)) if rl>0 else rl;print(tag,'rotated',n,float(normalized),flush=True)
   if rl>0:best=n;bestl=normalized
  if fixed:assert str(bestl)==str(old['lower']) and best==int(old['best']);print('fixed PASS',best, float(bestl))
  else:np.savez_compressed(f,C=C,R=R,scale=S,best=best,lower=str(bestl),allow=str(allow),L=np.array(L,str))
