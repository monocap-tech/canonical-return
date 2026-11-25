from bv_joint_trial import *
from itertools import product
# Direct response only, C=I. All remainder inequalities are Loewner inequalities.
def build(d,a):
 G=d['G'];Ge=rows(d['Ge']);Z=d['Z'];Ze=d['Ze']; assert d['ZS']==S*S;assert np.array_equal(a['ids'],d['c']['ids']);U=obj(a['U']);Ue=obj(a['Ue']);Gu=cf(obj(a['G']),S**3);Gue=rows(cf(obj(a['Ge']),S**3));T=obj(a['T']);ix=np.searchsorted(a['ids'],a['j']);Gvu=cf(Z[ix].T@T,S**3);Gvue=norm(Ze[ix].T@abs(T),S**3)
 Z2u=cf(obj(a['Z2']),S**4);Z2ue=rows(cf(obj(a['Z2e']),S**4));ut=Gu-Z2u;ute=Gue+Z2ue
 w=np.array([scalarfloor((l+h)/2) for l,h in zip(d['bl'],d['bh'])],object);we=max(max(abs(F(int(v),S)-l),abs(F(int(v),S)-h)) for v,l,h in zip(w,d['bl'],d['bh']));mid=(d['bmin']+d['bmax'])/2;rad=(d['bmax']-d['bmin'])/2;M=max(abs(d['bmin']),d['bmax']);b2max=M*M
 c,ce=cross(Z,Ze,U,Ue,w);t,te=cross(Z,Ze,U,Ue,np.ones(len(U),object));Vut=Gvu-cf(t,S**4);Vute=Gvue+norm(te,S**4);Bvu=cf(c,S**5)+mid*Vut;Bvue=norm(ce,S**5)+we*sqrtup((rows(d['Z2'])+rows(d['Z2E']))*(rows(Z2u)+Z2ue))+abs(mid)*Vute
 c,ce=gram(U,Ue,w);Bu=cf(c,S**5);Bue=rows(cf(ce,S**5))+we*(rows(Z2u)+Z2ue)
 w2=np.array([scalarceil(max(l*l,h*h)) for l,h in zip(d['bl'],d['bh'])],object);j2,j2e=gram(U,Ue,w2);J2=cf(j2,S**5);J2e=rows(cf(j2e,S**5))
 HG=cf(obj(a['HG']),S**3);HGe=rows(cf(obj(a['HGE']),S**3));CX=cf(obj(a['CX']),S**3);CE=rows(cf(obj(a['CE'])+obj(a['CE']).T,S**3));Y2=cf(obj(a['Y2']),S**4);Y2e=rows(cf(obj(a['Y2e']),S**4));ht=HG-Y2;hte=HGe+Y2e;O=CX+CX.T-HG;Oe=CE+HGe
 n=G.shape[0];eye=I(n);vt=d['Vtail'];vte=d['et'];base=G-Bvu-Bvu.T+Bu+J2-O+d['diag'];err=Ge+2*Bvue+Bue+J2e+Oe
 # PSD upper envelopes for actual tail Grams. The sign of the tail cross is handled by Young's inequality.
 V=vt+vte*eye;U=ut+ute*eye;H=ht+hte*eye
 base=base+err*eye+(d['bmax']+b2max)*U
 dd=d['derr']+d['e'];E=d['E']-d['eE']*eye
 print('base uniform finite error',float(err),'dd',float(dd),'tail errors',*[float(t) for t in (vte,ute,hte)],flush=True)
 return locals()
def compare(b,theta,kappa,gamma):
 H=b['base']+b['rad']*(theta*b['V']+b['U']/theta)+b['M']*(kappa*b['U']+b['H']/kappa)
 L=b['E']-((1+gamma)*H+(1+1/gamma)*b['dd']**2*b['eye'])/b['d']['mu'];return L,H
if __name__=='__main__':
 d=loadjoint();a=np.load(parent/'cr_int_t1d_ca_odd'/'odd_joint_Uaction.npz' if '--baseline' in sys.argv else root/'odd_joint_Uaction.npz');b=build(d,a);fixed='--fixed' in sys.argv; f=root/('cb_matrix_baseline.npz' if '--baseline' in sys.argv else 'cb_matrix_accepted.npz')
 if fixed:
  o=np.load(f);p=[F(str(v)) for v in o['parameters']]
 else:
  # numerical parameter search only. Use actual matrix Loewner upper bound, not uniform tails.
  E=np.array(b['E'],float);B=np.array(b['base'],float);V=np.array(b['V'],float);U=np.array(b['U'],float);H=np.array(b['H'],float);rad=float(b['rad']);M=float(b['M']);mu=float(d['mu']);dd=float(b['dd']);eye=np.eye(len(E));best=(-1e9,None)
  for th,ka,ga in product([F(1,j) for j in (10,20,30,40,50,60,80,100,150,200)],[F(j,2) for j in (1,2,3,4,6,8,12,16,24)],[F(1,j) for j in (50,80,100,125,150,175,200,250,300)]):
   hh=B+rad*(float(th)*V+U/float(th))+M*(float(ka)*U+H/float(ka));ll=E-((1+float(ga))*hh+(1+1/float(ga))*dd*dd*eye)/mu;v=np.linalg.eigvalsh(ll)[0]
   if v>best[0]:best=(v,[th,ka,ga])
  print('best proposal',best,flush=True);p=best[1]
 L,H=compare(b,*p);ql=rows(cf(d['q'].T@d['q'],S*S)-I(len(L)));low=exact_lower(L);nl=low/(1+ql) if low>0 else low;print('parameters',p,'head upper',float(rows(H)),'exact lower',float(nl),flush=True)
 if fixed:assert str(nl)==str(o['lower']) and nl>0;print('fixed PASS')
 else:np.savez_compressed(f,parameters=np.array([str(x) for x in p]),L=np.array(L,str),H=np.array(H,str),lower=str(nl),scale=S);print('PASS' if nl>0 else 'FAIL')
