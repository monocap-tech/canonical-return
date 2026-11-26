from pathlib import Path
p=Path(__file__).resolve().parent
s=(p.parent/'cr_alp_054/generate.py').read_text()
s=s.replace('for bl in bank[\'blocks\']:', 'for bl in bank[\'blocks\']:\n if len(sys.argv)>1 and bl[\'parity\']!=int(sys.argv[1]):continue')
pos=s.index(' # Complete finite trial space G plus outside pole-null D')
s=s[:pos]+''' # Append the genuine global bridge at cutoff 850, relative to Q600.
 # h = (c_D - r2 R600 b)/norm2; a finite sine combination plus complete b.
 sD=cc; sT=si; r2=sD/sT; mass2=sD*(sT-sD)/sT; norm2=mass2.sqrt()
 assert mass2.l>0 and (sT-sD).l>0
 fullids=allids+high
 d2=mat([[r2*c(j)/norm2] for j in allids]+[[c(j)/norm2] for j in high]); db2=-r2/norm2
 Bh=mul(d2.T,np.concatenate((mul(hlow,V),U),axis=0))+mul(hbp,V)*db2
 # Prime interactions against every sine and the complete pole.
 Pfull=mat([[pss(i,j) for j in fullids] for i in fullids]); pbfull=mat([[pb(j)] for j in fullids]); ph=mul(Pfull,d2)+pbfull*db2
 pbh=mul(pbfull.T,d2)[0,0]+Pbb*db2
 pgh=mul(d.T,ph[:len(allids),:])[0,0]+db*pbh
 phh=mul(d2.T,ph)[0,0]+db2*pbh
 # All gamma cross terms vanish by exact orthogonality to the older trial.
 rbridge=Bh-mul(ph[len(low):len(allids),:].T,X[:nt,:])-X[nt:,:]*pgh
 Gbridge=sym(mul(rbridge.T,rbridge))
 print('parity',k,'new bridge squared residual',float(sum(Gbridge[i,i].mid() for i in range(n))),flush=True)
 # Reuse the pinned 054 bank, which exactly represents the old J.
 raw054=(ROOT.parent/'cr_alp_054'/f'bank_{k}.json').read_bytes(); b54=json.loads(raw054)
 assert hashlib.sha256(raw054).hexdigest()==PIN54[k]
 RP=load(b54['R']); PP=load(b54['P']); KP=load(b54['K']);
 assert all(max(x.l,y.l)<=min(x.h,y.h) for x,y in zip(KP.flat,Knew.flat))
 assert all(max(x.l,y.l)<=min(x.h,y.h) for x,y in zip(RP.flat,np.concatenate((fa,Yp,Raug[nt:,:]),axis=0).flat))
 pa=ph[len(low):len(allids),:]; pd=ph[len(allids):,:]
 pa=pa-mul(ct,mat([[x/sa for x in mul(ct.T,pa)[0]]]))
 pd=pd-mul(ch,mat([[x/cc for x in mul(ch.T,pd)[0]]]))
 pcol=np.concatenate((pa,pd,mat([[pgh]])),axis=0)
 PP=np.block([[PP,pcol],[pcol.T,mat([[phh]])]])
 RP=np.concatenate((RP,rbridge),axis=0)
 out=dict(parity=k,scale=str(S),trial=600,band_end=N,predecessor_sha256=hashlib.sha256(raw054).hexdigest(),K=save(Knew),R=save(RP),P=save(PP),Gbridge=save(Gbridge),rbridge=save(rbridge),witness=saved['comparison_witness'])
 (ROOT/f'bank_{k}.json').write_text(json.dumps(out)+'\\n')
 print('parity',k,'complete E850 minus E85 residual bank',len(RP),'redundant rows',flush=True)
'''
import hashlib
pins={k:hashlib.sha256((p.parent/'cr_alp_054'/f'bank_{k}.json').read_bytes()).hexdigest() for k in (1,2)}
s=s.replace('S=10**100;',f'PIN54={pins!r}\nS=10**100;')
(p/'generate.py').write_text(s)
