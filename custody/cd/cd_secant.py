from pathlib import Path
import sys, hashlib, json
from fractions import Fraction as F
workspace=Path(__file__).resolve().parent
par=sys.argv[1]; assert par in ('even','odd')
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
sys.path.insert(0,str(source))
from bv_joint_trial import *
code=(workspace/'cr_int_t1d_cb'/'cb_matrix.py').read_text();exec(code[:code.index("if __name__=='__main__':")])
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
for path,manifest in [(source,'ca_sha256.txt' if par=='even' else 'cb_sha256.txt'),(workspace/'cr_int_t1d_cb','cb_sha256.txt'),(workspace/'cr_int_t1d_bp','bp_sha256.txt'),(workspace/'cr_int_t1d_bs','bs_sha256.txt'),(workspace/'cr_int_t1d_bt','bt_sha256.txt'),(workspace/'cr_int_t1d_bu','bu_sha256.txt')]:
 for line in (path/manifest).read_text().splitlines():
  h,name=line.split(maxsplit=1);assert hashlib.file_digest((path/name.strip()).open('rb'),'sha256').hexdigest()==h,(path,name)
d=loadjoint();a=np.load(source/f'{par}_joint_Uaction.npz');b=build(d,a);n=len(b['E']);eye=b['eye'];Z=d['Z'];Ze=d['Ze'];U=obj(a['U']);Ue=obj(a['Ue']);Y=obj(a['Y']);Ye=obj(a['YE'])
m=b['mid'];r=b['rad'];s=m+F(1,2);c0=m*m+m+r*r
assert c0>=0 and r>=0
# Actual full u*h and v*h vanish. The finite Fourier synthesis centers below
# are rational, so divide exactly rather than casting fractional entries to int.
ww,wwe=gram(Z+Y-s*U,Ze+Ye+abs(s)*Ue)
Gu=b['Gu'];Gvu=b['Gvu'];full=d['G']+b['HG']+s*s*Gu-s*(Gvu+Gvu.T)
err=rows(d['Ge'])+b['HGe']+s*s*b['Gue']+2*abs(s)*b['Gvue']+rows(wwe/S**4)
W=full-ww/S**4+err*eye
assert np.array_equal(W,W.T)
c,ce=cross(U,Ue,Y,Ye,np.ones(len(U),object));midc=cf(c,S**4);mide=norm(ce,S**4)
base=b['base']-(d['bmax']+b['b2max'])*b['U']+c0*b['U']+m*(midc+midc.T)+2*abs(m)*mide*eye
Q=cf(d['q'].T@d['q'],S*S);X=cf(d['x'].T@d['x'],S**4);aq=d['eY']+sqrtup(d['tau']);ax=d['eU']+sqrtup(d['tauU']);er=d['e'];nq=sqrtup(rows(Q));nx=sqrtup(rows(X));v=[aq*nq,ax*nx,er];total=sum(v);lamb=[vi/total for vi in v];assert all(li>0 for li in lamb) and sum(lamb)==1
D2=aq**2/lamb[0]*Q+ax**2/lamb[1]*X+er**2/lamb[2]*eye
eta=nx/nq;assert eta>0
Edir=d['E']-d['eK']*(eta*Q+X/eta)-d['eD']*X
ql=rows(Q-eye)
# Search on a frozen rational grid by float only; exact LDL supplies final acceptance.
if '--fixed' in sys.argv:
 old=json.loads((workspace/'cd_results'/f'{par}.json').read_text());th,ga=map(F,old['parameters'])
else:
 best=(-float('inf'),None);Ef=np.array(Edir,float);Bf=np.array(base,float);Wf=np.array(W,float);Uf=np.array(b['U'],float);Df=np.array(D2,float);rf=float(r);mu=float(d['mu'])
 for ti in range(5,121):
  theta=F(1,ti);H=Bf+rf*(float(theta)*Wf+Uf/float(theta))
  for gi in range(40,321,5):
   gamma=F(1,gi);L=Ef-((1+float(gamma))*H+(1+1/float(gamma))*Df)/mu;ev=np.linalg.eigvalsh(L)[0]
   if ev>best[0]:best=(ev,(theta,gamma))
 th,ga=best[1];print('float proposal',best,flush=True)
H=base+r*(th*W+b['U']/th);L=Edir-((1+ga)*H+(1+1/ga)*D2)/d['mu'];low=exact_lower(L)/(1+ql);assert low>0
margin=F(int(low*100000//1),100000);assert margin<low
A=L-margin*(1+ql)*eye;assert np.array_equal(A,A.T);T=np.vectorize(lambda v:int(v*S//1),otypes=[object])(A)-n*eye;assert ldl_positive(cf(T,S))
result=dict(parity=par,n=n,parameters=[str(th),str(ga)],lower=str(low),margin=str(margin),qdefect=str(ql),upper_triangle=[int(T[i,j]) for i in range(n) for j in range(i,n)])
out=workspace/'cd_results';out.mkdir(exist_ok=True)
if '--fixed' in sys.argv:assert result==old
else:(out/f'{par}.json').write_text(json.dumps(result,indent=2)+'\n')
print('exact PASS',par,n,str(low),float(low),'theta,gamma',th,ga,'W row upper',float(rows(W)),flush=True)
