from pathlib import Path
import sys,hashlib,json
from fractions import Fraction as F
workspace=Path(__file__).resolve().parent
par=sys.argv[1];assert par in ('even','odd')
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
sys.path.insert(0,str(source))
from bv_joint_trial import *
# use the archived CB build for both parities, with the actual parity-specific loader above.
from importlib.machinery import SourceFileLoader
# avoid importing CB's odd loader in the even audit; execute only common functions.
code=(workspace/'cr_int_t1d_cb'/'cb_matrix.py').read_text();exec(code[:code.index("if __name__=='__main__':")])
source=workspace/('cr_int_t1d_ca' if par=='even' else 'cr_int_t1d_cb')
for path,manifest in [(source,'ca_sha256.txt' if par=='even' else 'cb_sha256.txt'),(workspace/'cr_int_t1d_cb','cb_sha256.txt'),(workspace/'cr_int_t1d_bp','bp_sha256.txt'),(workspace/'cr_int_t1d_bs','bs_sha256.txt'),(workspace/'cr_int_t1d_bt','bt_sha256.txt'),(workspace/'cr_int_t1d_bu','bu_sha256.txt')]:
 for line in (path/manifest).read_text().splitlines():
  h,name=line.split(maxsplit=1);assert hashlib.file_digest((path/name.strip()).open('rb'),'sha256').hexdigest()==h,(path,name)
d=loadjoint();a=np.load(source/f'{par}_joint_Uaction.npz');b=build(d,a);n=len(b['E']);eye=b['eye'];Z=d['Z'];Ze=d['Ze'];U=obj(a['U']);Ue=obj(a['Ue']);Y=obj(a['Y']);Ye=obj(a['YE'])
# Full physical v and u have support in I; h has support in the exterior. Their full cross Grams vanish.
ww,wwe=gram(Z+Y,Ze+Ye);W=d['G']+b['HG']-cf(ww,S**4)+(rows(d['Ge'])+b['HGe']+rows(cf(wwe,S**4)))*eye
c,ce=cross(U,Ue,Y,Ye,np.ones(len(U),object));midc=cf(c,S**4);mide=norm(ce,S**4);base=b['base']+b['mid']*(midc+midc.T)+2*abs(b['mid'])*mide*eye
th,ga=(F(1,25),F(1,80)) if par=='even' else (F(1,50),F(1,175))
H=base+b['rad']*(th*W+b['U']/th);L=b['E']-((1+ga)*H+(1+1/ga)*b['dd']**2*eye)/d['mu'];ql=rows(cf(d['q'].T@d['q'],S*S)-I(n));low=exact_lower(L);nl=low/(1+ql);assert nl>0
expected=F('51670000000000000000000/166666666666972004399661') if par=='even' else F('46070000000000000000000/1000000000006244646926501');assert nl==expected
# Directed, directional error upper bounds; weights via the already certified scalar norms.
Q=cf(d['q'].T@d['q'],S*S);X=cf(d['x'].T@d['x'],S**4);aq=d['eY']+sqrtup(d['tau']);ax=d['eU']+sqrtup(d['tauU']);er=d['e'];nq=sqrtup(rows(Q));nx=sqrtup(rows(X));v=[aq*nq,ax*nx,er];total=sum(v);lamb=[vi/total for vi in v];assert all(li>0 for li in lamb) and sum(lamb)==1
D2=aq**2/lamb[0]*Q+ax**2/lamb[1]*X+er**2/lamb[2]*eye
eta=nx/nq;assert eta>0
Edir=d['E']-d['eK']*(eta*Q+X/eta)-d['eD']*X
Ld=Edir-((1+ga)*H+(1+1/ga)*D2)/d['mu'];ld=exact_lower(Ld)/(1+ql)
# Preserve conservative shifted integer comparator for independent exact verification of each strict margin.
result=[]
for name,M,margin in [('joint',L,F(31,100) if par=='even' else F(46,1000)),('directional',Ld,F(int(ld*100000//1),100000))]:
 assert margin>0 and margin<(nl if name=='joint' else ld)
 A=M-margin*(1+ql)*eye;assert np.array_equal(A,A.T)
 T=np.vectorize(lambda v:int(v*S//1),otypes=[object])(A)-n*eye;assert ldl_positive(cf(T,S))
 result.append(dict(name=name,parity=par,n=n,margin=str(margin),qdefect=str(ql),upper_triangle=[int(T[i,j]) for i in range(n) for j in range(i,n)]))
out=workspace/'cc_results';out.mkdir(exist_ok=True);(out/f'{par}_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print('fixed joint PASS',par,n,str(nl),float(nl),'directional',str(ld),float(ld),flush=True)
print('joint W bound',float(rows(W)),'midpoint cross error',float(mide),'directional weights',*[float(x) for x in lamb],flush=True)
