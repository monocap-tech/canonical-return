from bv_response import *
from math import factorial
# exact directed Gram primitives from BU/BS
src=(parent/'cr_int_t1d_bs'/'bs_circle.py').read_text();r=isqrt(2*S*S);exec(src[src.index('atan=lambda'):src.index('# g_j')]);src=(parent/'cr_int_t1d_bu'/'bu_circle.py').read_text();exec(src[src.index('qlo='):src.index('# directed cosh')]);even_gram=gram
# restore generic weighted gram
from bv_response import gram as weighted_gram
def odd_gram(j,k):
 d=j[:,None]-k[None,:];s=j[:,None]+k[None,:];l,h=term(d);u,v=term(s);return l-v+(d==0)*(S//2),h-u+(d==0)*(S//2)
if __name__=='__main__':
 par=sys.argv[1];tag='surplus' if par=='odd' else 'raw';g=odd_gram if par=='odd' else even_gram;d=load(par,tag);a=np.load(root/f'{par}_joint_proposal.npz');j=a['j'];T=a['T'];ids=d['c']['ids'];ix=np.searchsorted(ids,j);assert np.array_equal(ids[ix],j);UU=[];UE=[]
 for start in range(0,len(ids),256):
  l,h=g(ids[start:start+256],j);mid=(l+h)//2;rad=np.maximum(mid-l,h-mid);assert rad.max()<=5;UU.append(dot(mid,T));UE.append(dot(rad,abs(T)))
 U=np.concatenate(UU);Ue=np.concatenate(UE);G=T.astype(object).T@U[ix];Ge=abs(T).astype(object).T@Ue[ix];G=(G+G.T)//2;Ge=np.maximum(Ge,Ge.T)+1
 # actual circle J*Eu enclosure at common S
 bl=np.array([scalarfloor(l) for l in d['bl']],object);bh=np.array([scalarceil(h) for h in d['bh']],object);al=np.minimum.reduce([b[:,None]*v for b in [bl,bh] for v in [U-Ue,U+Ue]])//S**2;ah=-((-np.maximum.reduce([b[:,None]*v for b in [bl,bh] for v in [U-Ue,U+Ue]]))//S**2);A=((al+ah)//2).astype(np.int64);AE=np.maximum(A.astype(object)-al,ah-A.astype(object));assert np.max(AE)<100000
 # exterior trial P_out of largest action coefficients, exact upper/lower Gram
 n=int(sys.argv[2]);ind=np.sort(np.argsort(np.sum(A.astype(float)**2,axis=1))[-n:]);k=ids[ind];V=A[ind];HH=[];HE=[]
 for start in range(0,len(ids),256):
  jj=ids[start:start+256];l,h=g(jj,k);delta=(jj[:,None]==k[None,:])*S;l,h=delta-h,delta-l;mid=(l+h)//2;rad=np.maximum(mid-l,h-mid);assert rad.max()<=5;HH.append(dot(mid,V));HE.append(dot(rad,abs(V)))
 Y=np.concatenate(HH);YE=np.concatenate(HE);HG=V.astype(object).T@Y[ind];HGE=abs(V).astype(object).T@YE[ind];HG=(HG+HG.T)//2;HGE=np.maximum(HGE,HGE.T)+1;CX=A.astype(object).T@Y;CE=abs(A).astype(object).T@YE+AE.T@abs(Y)+AE.T@YE
 Z2,Z2e=weighted_gram(U,Ue);Y2,Y2e=weighted_gram(Y,YE)
 np.savez_compressed(root/f'{par}_joint_Uaction.npz',ids=ids,j=j,T=T,U=np.array(U,str),Ue=np.array(Ue,str),G=np.array(G,str),Ge=np.array(Ge,str),A=A,AE=np.array(AE,str),idx=ind,V=V,Y=np.array(Y,str),YE=np.array(YE,str),HG=np.array(HG,str),HGE=np.array(HGE,str),CX=np.array(CX,str),CE=np.array(CE,str),Z2=np.array(Z2,str),Z2e=np.array(Z2e,str),Y2=np.array(Y2,str),Y2e=np.array(Y2e,str),scale=S)
 print('exact directed Uaction PASS',par,len(ids),len(j),n,'max A error',np.max(AE),flush=True)
