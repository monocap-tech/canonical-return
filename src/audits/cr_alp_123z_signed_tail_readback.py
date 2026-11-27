#!/usr/bin/env python3
"""Signed pole-tail readback and independent sparse high-mode checks.
Frozen E850 correction and moments; no new joint solve.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_123y_signed_pole_tail import leading
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def digest(b):return hashlib.sha256(b).hexdigest()
def up(v):return max(abs(F(x)) for x in v)
def root(v):
    assert v>=0
    D=10**110;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def namespace():
    src=Path(__file__).with_name('cr_alp_046_full_finite.py')
    raw=src.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()=='554e7ab54da4545987d89aa660a4abe6c951cb63'
    s=raw.decode()
    def replace(a,b):
        nonlocal s
        assert s.count(a)==1,a
        s=s.replace(a,b)
    replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=2')
    replace('def psi_pair(x,y=ZERO):','''def atan_wide(x):
 z=ZERO
 for _ in range(8):
  if x.up()<=F(1,3):return z+atan(x)
  x0=x;x=(x-F(1,4))/(1+x/4);z+=atan(I(F(1,4)))
  assert x.l>=0 and x.up()<x0.up()
 raise AssertionError("atan range reduction")
def psi_pair(x,y=ZERO):''')
    replace('assert y.up()/F(xx.l,S)<F(1,3)','assert y.up()/F(xx.l,S)<5')
    replace('li=atan(y/xx)','li=atan_wide(y/xx)')
    replace('range(1,N+1)};trig=', '(1,2,2049,2050,2111,2112,3071,3072)};trig=')
    replace('for j in ids};Hb={}\n for j in ids:', 'for j in nu if j%2==k%2};Hb={}\n for j in c:')
    replace('return order,[[convert(H,i,j) for j in order],','return order,[[convert(H,i,j) for j in order],') if False else None
    replace(' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]',
            ' return c,Hb')
    g={'__file__':str(src),'__name__':'cr_alp_123z_sparse'}
    exec(compile(s,str(src),'exec'),g)
    return g

def check_poles(g,k,J,fb,T,H0,receipt):
    I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    def iv(v):return I(*map(F,v))
    A=iv(receipt['constant']);B=iv(receipt['log_coefficient'])
    C=list(map(iv,receipt['cosine_coefficients']))
    p2,p3=F(receipt['pointwise_inverse_square_upper']),F(receipt['pointwise_inverse_cube_upper'])
    rt,eh,em,qt=(g[x] for x in ('rt','eh','em','qt'))
    b0=(eh+em)/2 if k==1 else -(eh-em)/2
    b1=-(eh-em)/4 if k==1 else (eh+em)/4
    ep,en=(b0+2*b1)/2,(b0-2*b1)/2
    c,Hb=g['build'](k);rows=[]
    for j in c:
        if j<=J:continue
        v=g['nu'][j]
        P=ZERO
        for n,w,(si,co) in zip(g['NS'],g['W'],g['trig'][j]):
            sn=I(n).sqrt();bt=ep*sn+en/sn
            P-=2*w/(rt*(v.sq()+qt))*(v*(bt+b0*co)-b1*si)
        exact=H0*v/(v.sq()+qt)+fb*(Hb[j]-P)+T*P
        approx=(A+B*g['log'](I(F(j,J-1)))+sum((z*tr[1] for z,tr in zip(C,g['trig'][j])),ZERO))/v
        err=(exact-approx).up();vmin=F(v.l,S)
        bound=p2/vmin**2+p3/vmin**3
        assert err<=bound,(k,j,float(err),float(bound))
        rows.append({'index':j,'error_upper':str(err),'bound':str(bound)})
        print('pointwise',k,j,'ratio',float(err/bound) if bound else 0,flush=True)
    assert len(rows)==3
    return rows

def main(second,joint,quant,out):
    raw=second.read_bytes();assert digest(raw)=='9152dcb759260d498c7496e579285d83e28172b6bcb3192c32e961886335253f'
    d=json.loads(raw);assert d['status']=='CONDITIONAL_SECOND_AUXILIARY_CORRECTION'
    jraw=joint.read_bytes();assert digest(jraw)=='09b34469b32c13ebf5f5b690d1eb9ad3a9fb5676f171a1c7f85325cedc199d7b'
    g=json.loads(jraw)
    qraw=quant.read_bytes();assert digest(qraw)=='f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de'
    q=json.loads(qraw)
    ng=namespace();I=ng['I']
    J=2048;V=F(15*J,11);dc=F('1.822938456')
    for k,rows in d['sectors'].items():
        assert len(rows)==1
        r=rows[0];assert r['cutoff']==850 and r['source_cutoff']==365 and r['finite_prefix']==J
        w=next(x['joint_trial'] for x in g['sectors'][k]['rows'] if x['cutoff']==365)
        ds={int(i):up(v) for i,v in w['source_sines'].items()}
        fb=up(r['source_pole']);T=up(r['combined_pole'])
        # tau = fb - T, signed endpoints.
        tl,th=map(F,r['source_pole']);ul,uh=map(F,r['combined_pole'])
        tau=max(abs(tl-uh),abs(th-ul))
        ts={i:v+tau*F(33,15*i) for i,v in ds.items()}
        ts.update({int(i):abs(F(v)) for i,v in r['corrected_tail_coordinates'].items()})
        assert max(ts)<=850 and all(i>85 for i in map(int,r['corrected_tail_coordinates']))
        n=32;rho=(F(10*max(ts),7)/V)**2;assert rho<1
        coeff={}
        def charge(s,v):
            assert v>=0
            coeff[s]=coeff.get(s,F(0))+v
        def series(label,mult,shift,H):
            mom=r['signed_moments'][label];assert len(mom)==n
            for l,z in enumerate(mom):charge(2*l+2-shift,mult*up(z))
            charge(2*n+2-shift,mult*H/(1-rho))
        # |AG|,|AP|<4; c_j <= 33/(15j). Use independent rational
        # upper absolute moments, not the lost original remainder coefficients.
        H=sum(4*(ds.get(i,F(0))+t)*(F(10*i,7)**2+F(1,4))*F(10*i,7)**(2*n) for i,t in ts.items())
        series('centered_combined',F(20,11),-1,H)
        series('gamma_sine',F(80,11),0,sum(v*F(10*i,7)**(2*n+1) for i,v in ds.items()))
        series('prime_sine',F(80,11),0,sum(v*F(10*i,7)**(2*n+1) for i,v in ts.items()))
        # Count the actual parity in each sine-moment norm.
        J0=J-1
        series_tail=sum(z*root(F(J0,2*(2*s-1)))/F(15*J0,11)**s for s,z in coeff.items())
        pole_args=[I(*map(F,r[name])) for name in ('source_pole','combined_pole','leading_combined_moment')]
        pole,receipt=leading(ng,int(k),J,*pole_args)
        receipt['sparse_pointwise_checks']=check_poles(ng,int(k),J,*pole_args,receipt)
        receipt['sine_series_norm_upper']=str(series_tail)
        r['signed_pole_tail_receipt']=receipt
        tail=series_tail+pole
        U=F(r['weighted_prefix_upper'])+tail*tail/dc # omit pole subtraction safely
        E=F(r['joint_energy_lower']);b=F(r['finite_corrected_residual_norm_upper'])
        pp=q['sectors'][k];delta=F(1,10**5) if k=='1' else F(21,10**5)
        A=1/F(pp['lambda_K_lower'])+1/delta;B=F(pp['kappa_upper'])/delta
        loss=U+(root(A)*b+root(B*U))**2
        old=F(r['infinite_tail_norm_upper'])
        r.update({'unrepaired_infinite_tail_norm_upper':str(old),'infinite_tail_norm_upper':str(tail),
                  'weighted_inverse_error_upper':str(U),'weighted_pole_subtraction_lower':'0',
                  'response_loss_upper':str(loss),'conditional_pole_schur_lower':str(E-loss),
                  'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'})
        print(k,'signed tail',float(tail),'prefix U',float(F(r['weighted_prefix_upper'])),
              'complete U',float(U),'U/E',float(U/E) if E>0 else None,
              'loss/E',float(loss/E) if E>0 else None,flush=True)
    d['repair']={'original_sha256':digest(raw),'method':'signed pole asymptotic and parity Abel estimate, independent sparse checks, parity sine series integral; omit nonnegative weighted pole subtraction'}
    out.write_text(json.dumps(d,indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
