from pathlib import Path
p=Path(__file__).resolve().parent;s=(p.parent/'cr_alp_055/generate.py').read_text(); head=s.split('nu={j:')[0].replace('N=850;KTRIAL=600','N=1050;KTRIAL=600').replace('assert y.up()/F(xx.l,S)<F(3,5)','assert y.up()/F(xx.l,S)<F(4,5)').replace('li=(atan(y/xx) if (y/xx).up()<=F(1,3) else atan(I(F(1,4)))+atan((y/xx-F(1,4))/(1+y/xx/4)))','li=atan_wide(y/xx)');head=head.replace('def psi_pair(x,y=ZERO):','def atan_wide(t):\n z=ZERO\n while t.up()>F(1,3):\n  z+=atan(I(F(1,4)));t=(t-F(1,4))/(1+t/4)\n return z+atan(t)\ndef psi_pair(x,y=ZERO):')
part=s[s.index('def hss(i,j):'):s.index(' band=list(range(87')]; part=part.replace('for j in nu}', 'for j in nu}')
loop='''nu={j:j*NUF for j in range(1,N+1)};trig={};R={};G={};AI={};R2={}
for j,v in nu.items():
 trig[j]=[sincos(v*t) for t in TS]
 if j<=85 or j>850:
  re,im,tre=psi_pair(qt,v/2);R[j]=im/(2*v);G[j]=re-log(PI);AI[j]=v*(R[j]-(-1)**j*sum(e/(b+v.sq()) for e,b in zip(EX,B2)).widen(1))+sum(w*t[0] for w,t in zip(W,trig[j]))
 if j%100==0:print('coefficients',j,flush=True)
'''
end='''
 # Frozen trial in F600 plus its bridge. New correction in F600, D850, g, h.
 low=list(range(k,86,2)); Aids=list(range(87 if k==1 else 86,601,2)); Dids=list(range(601 if k==1 else 602,851,2)); fullids=low+Aids+Dids; na=len(Aids); nd=len(Dids)
 def col(js):return mat([[c(j)] for j in js])
 cl,ca,cd=map(col,(low,Aids,Dids));sl=mul(cl.T,cl)[0,0];sa=mul(ca.T,ca)[0,0];sd=mul(cd.T,cd)[0,0];st=BM-sl;si=st-sa;r=sa/st;r2=sd/si;norm=(sa*(st-sa)/st).sqrt();norm2=(sd*(si-sd)/si).sqrt();db=-r/norm;db2=-r2/norm2
 d=mat([[r*c(j)/norm] for j in low]+[[c(j)/norm] for j in Aids]+[[ZERO] for j in Dids]);d2=mat([[r2*c(j)/norm2] for j in low+Aids]+[[c(j)/norm2] for j in Dids]);
 pin=(ROOT.parent/'cr_alp_053'/f'result_{k}.json').read_bytes();assert hashlib.sha256(pin).hexdigest()=={1:'a133a5bf9a8229b16e6459aee952bdca677567ed9df425e7430c466bd87e0877',2:'ee7ae0d12d7beefd2ecb08c3825ca5f7ab9ea92da2fba2f28b836094c212c4f3'}[k];saved=json.loads(pin);free=load(saved['trial_free_rows']);last=-mul(mat([[c(j)/c(Aids[-1]) for j in Aids[:-1]]]),free);X=np.concatenate((free,last),axis=0);xg=load(saved['bridge_coeff'])
 cor=json.loads((ROOT/f'correction_{k}.json').read_text()); sc=int(cor['scale']);assert S%sc==0;Z=mat([[I(int(l)*(S//sc),int(h)*(S//sc),True) for l,h in row] for row in cor['Z']]);assert len(Z)==na+nd+2
 za=Z[:na,:];zd=Z[na:na+nd,:];za=za-mul(ca,mat([[v/sa for v in mul(ca.T,za)[0]]]));zd=zd-mul(cd,mat([[v/sd for v in mul(cd.T,zd)[0]]]));zg=Z[-2:-1,:]+xg;zh=Z[-1:,:]
 D=np.concatenate((mat([[ZERO]*n for j in low]),X+za,zd),axis=0)+mul(d,zg)+mul(d2,zh);Db=zg*db+zh*db2
 high=list(range(851 if k==1 else 852,N+1,2)); C=col(high)
 U=mul(mat([[hss(i,j)-c(i)/c(k)*hss(k,j) if i!=k else hss(k,j)-c(k)/BM*hb(j) for i in ids] for j in high]),V)
 Pout=mat([[pss(j,i) for i in fullids] for j in high]);Y=U-mul(Pout,D)-mul(mat([[pb(j)] for j in high]),Db)
 # gamma I contribution vanishes after Q850, since X+Z lies in E850.
 cc=mul(C.T,C)[0,0];Yp=Y-mul(C,mat([[v/cc for v in mul(C.T,Y)[0]]])); Gout=sym(mul(Yp.T,Yp))
 K=mat([[I(int(l)*(S//sc),int(h)*(S//sc),True) for l,h in row] for row in cor['K']]);GI=mat([[I(int(l)*(S//sc),int(h)*(S//sc),True) for l,h in row] for row in cor['G']]);
 Cmp=sym(K-(GI+Gout)*I(F(200,51)))
 def wit(M,largest=False):
  _,vs=np.linalg.eigh(midfloat(M));z=[round(float(v)*10**12) for v in vs[:, -1 if largest else 0]]; zz=mat([[I(v)] for v in z]);val=mul(zz.T,mul(M,zz))[0,0]/sum(v*v for v in z);return {'z':z,'interval':compact(val),'lower':str(F(val.l,S)),'upper':str(F(val.h,S))}
 nw=wit(Gout,True);cw=wit(Cmp);print(k,'outside squared norm lower',float(F(nw['lower'])),'comparison upper',float(F(cw['upper'])),flush=True)
 out=dict(parity=k,scale=str(S),end=N,Gout=save(Gout),K=save(K),Gin=save(GI),norm_witness=nw,comparison_witness=cw,correction_sha256=hashlib.sha256((ROOT/f'correction_{k}.json').read_bytes()).hexdigest());(ROOT/f'result_{k}.json').write_text(json.dumps(out)+'\\n')
'''
(p/'audit.py').write_text(head+loop+part+end)
