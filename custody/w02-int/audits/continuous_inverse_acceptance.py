"""Independent arithmetic acceptance replay for W02 INT continuous inverse.
Run with numpy and sympy:
  python continuous_inverse_acceptance.py --identity
  python continuous_inverse_acceptance.py receiver.npz comparison.npz result.json
The second mode checks exact assembly and positivity, not the analytic/source
hypotheses or the underlying directed integral evaluation.
"""
import sys, json, hashlib
from fractions import Fraction as F
from pathlib import Path

def identity():
    import sympy as s
    for N in (4,8):
        j=[1,N+1,2*N+1,N-1,2*N-1,N,2*N,3*N+2]
        C=s.Matrix([[i%3-1,(i*i+1)%5-2] for i in range(len(j))])
        A=s.zeros(N,2); B=s.zeros(N,2); V=s.zeros(N,2)
        for l,k in enumerate(j):
            phase=-s.pi*k/2+s.pi*k/N
            sc=s.sin(s.pi*k/N)/(s.pi*k/N)
            for c in range(2):
                A[k%N,c]+=C[l,c]*s.cos(phase)*sc
                B[k%N,c]+=C[l,c]*s.sin(phase)*sc
            for m in range(N):
                x0=-4+s.Rational(16*m,N);x1=x0+s.Rational(16,N)
                w=s.pi*k/8
                for c in range(2):
                    V[m,c]+=C[l,c]*(s.cos(w*x0)-s.cos(w*x1))/(2*w*s.Rational(16,N))
        X=s.Matrix(N,2,lambda k,c:A[k,c]-A[-k%N,c])
        Y=s.Matrix(N,2,lambda k,c:B[k,c]+B[-k%N,c])
        assert all(s.simplify(z)==0 for z in X.T*X+Y.T*Y-s.Rational(16,N)*V.T*V)
        print('exact symbolic alias/reflection Parseval identity PASS N=',N)

def replay(receiver,comparison,result):
    import numpy as np
    r=json.loads(Path(result).read_text())
    assert r['status']=='PASS'
    assert r['receiver_evidence']=='a4b63576f4a617c42975fdae85898fd1842739d6'
    assert hashlib.sha256(Path(receiver).read_bytes()).hexdigest()=='834d894be2bdac4bf7851d18862b203e7b234a68b6116f7d7d00a9c781929dac'
    assert hashlib.sha256(Path(comparison).read_bytes()).hexdigest()==r['comparison']['output_sha256']
    assert r['fft_impulse_control']==r['aliased_cell_integral_control']=='PASS'
    assert r['continuous']['exact_shift_overlay']==r['continuous']['integer_fft']=='PASS'
    v=np.load(receiver,allow_pickle=False); a=np.load(comparison,allow_pickle=False)
    f=lambda x:np.vectorize(lambda y:F(str(y)),otypes=[object])(x)
    I=np.eye(32,dtype=object)
    d,t,mu=(F(str(v[k])) for k in ('delta','t','mu'))
    assert d>0 and t>0 and mu==F(53,10000)
    for k in ('Gu','Gz','Gs','Gr','Fup','H','L'):
        z=f(a[k]);assert z.shape==(32,32) and np.array_equal(z,z.T)
    H=f(a['Fup'])+(1+d*(1+1/t))*f(a['Gr'])+d*(1+t)*f(a['Gu'])
    assert np.array_equal(H,f(a['H']))
    assert F(str(v['energy_error']))>=0
    L=f(v['E'])-f(v['penalty'])-H/mu-F(str(v['energy_error']))*I
    assert np.array_equal(L,f(a['L']))
    S=int(a['scale'])
    A=np.vectorize(lambda x:int(x*S//1),otypes=[object])(L)-32*I
    assert np.array_equal(A,np.vectorize(int,otypes=[object])(a['integer_comparison']))
    qdef=F(str(a['Q_gram_defect']));assert 0<=qdef<1
    assert qdef==F(r['comparison']['Q_gram_defect'])
    lower=F(str(a['lower']));assert lower==F(r['comparison']['lower'])
    # Independent fraction-free Sylvester test of shifted leading principal
    # minors. For symmetric matrices strict positive leading minors imply SPD.
    B=np.vectorize(int,otypes=[object])(A)
    assert (lower*S).denominator==1
    B=B-int(lower*S)*I
    prev=1; positive=True
    for k in range(32):
        pivot=int(B[k,k])
        if pivot<=0:
            positive=False;break
        for i in range(k+1,32):
            for j in range(k+1,32):
                x=int(B[i,j])*pivot-int(B[i,k])*int(B[k,j])
                q,rem=divmod(x,prev);assert rem==0
                B[i,j]=q
        prev=pivot
    assert positive==r['comparison']['exact_ldl_positive']
    if positive: assert lower>=0
    else: assert lower==0
    print('exact assembly and independent Sylvester replay PASS',
          'strict positivity' if positive else 'insufficient bound',
          'lower',lower)
    print('Analytic and source hypotheses remain separate obligations.')

if __name__=='__main__':
    if sys.argv[1:]==['--identity']:identity()
    else:
        assert len(sys.argv)==4
        replay(*sys.argv[1:])
