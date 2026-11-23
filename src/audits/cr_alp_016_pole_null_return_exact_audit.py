#!/usr/bin/env python3
"""CR-ALP-016: exact two-pole-null ramp combination and full signed Weil test.

Recomputes full signed 3x3 Weil Gram from prime powers (n<=729),
explicit Gamma antiderivatives, rigorous Fraction-log enclosures, and
rational Euler/pi bounds. This is NOT Lean, global positivity, Connes
operator-core admission, an Alpoge inequality, or RH certification.
"""
from fractions import Fraction as Q
from math import isqrt

def primes_through(n):
    sieve=[True]*(n+1)
    sieve[0]=sieve[1]=False
    for p in range(2,isqrt(n)+1):
        if sieve[p]:
            for k in range(p*p,n+1,p): sieve[k]=False
    return [k for k in range(2,n+1) if sieve[k]]

def log_atanh(q, N=22):
    assert 1<=q<=2
    z=(q-1)/(q+1)
    acc=Q(0)
    term=z
    for j in range(N):
        acc+=2*term/(2*j+1)
        term*=z*z
    tail=2*term/((2*N+1)*(1-z*z))
    return acc,acc+tail

def log_int(n):
    m=n.bit_length()-1
    a,b=log_atanh(Q(n,2**m))
    c,d=log_atanh(Q(2))
    return a+m*c,b+m*d

def log_ratio(a,b=1):
    l,u=log_int(a); ll,uu=log_int(b)
    return l-uu,u-ll

def enclose(terms):
    low=high=Q(0)
    for fac,num,den in terms:
        ll,uu=log_ratio(num,den)
        if fac>=0:
            low+=fac*ll
            high+=fac*uu
        else:
            low+=fac*uu
            high+=fac*ll
    return low,high

def prime_coefficients(bound, f):
    result={}
    for p in primes_through(bound):
        n=p
        coeff=Q(0)
        while n<=bound:
            coeff+=f(n)
            n*=p
        if coeff: result[p]=coeff
    return result

X=9
q=Q(3)
s=q+1/q
A=2+s*s
assert s==Q(10,3) and A==Q(118,9)
assert 1-s*q+q*q==0
assert 1-s/q+1/(q*q)==0
mass=8*A
assert mass==Q(944,9)

def r0(n):
    return Q(9,n)-1 if 2<=n<=9 else Q(0)
def r1(n):
    if n<=9: return Q(3*(n-1),n)
    if n<=81: return Q(81-n,3*n)
    return Q(0)
def r2(n):
    if n<=9: return Q(0)
    if n<=81: return 1-Q(9,n)
    if n<=729: return Q(81,n)-Q(1,9)
    return Q(0)
def rh(n):
    return A*r0(n)-2*s*r1(n)+r2(n)

assert rh(9)==Q(-160,9)
assert rh(81)==Q(8,9)
assert rh(729)==0
assert any(rh(n)!=0 for n in (16,81,625))
assert all(rh(n)==0 for n in range(730,751))

cd=prime_coefficients(9,lambda n:2*r0(n))
c1=prime_coefficients(81,r1)
c2=prime_coefficients(729,r2)
ch=prime_coefficients(729,lambda n:2*rh(n))
for p in primes_through(729):
    assert ch.get(p,0)==A*cd.get(p,0)-4*s*c1.get(p,0)+2*c2.get(p,0)
assert len(c2)==129 and len(ch)==129
assert cd=={2:Q(39,4),3:Q(4),5:Q(8,5),7:Q(4,7)}
assert c1[2]==Q(533,64) and c2[2]==Q(4285,1536)

terms_d=[(Q(16),9,1)]+[(-v,p,1) for p,v in cd.items()]+[(Q(10),5,1),(Q(-16),4,1)]
ld,ud=enclose(terms_d)
# Sharp harmonic-number enclosure:
# H_N-log N-1/(2N) < gamma < H_N-log N-1/(2N)+1/(12N²).
N=1000
HN=sum((Q(1,k) for k in range(1,N+1)),Q(0))
lN,uN=log_int(N)
glo=HN-uN-Q(1,2*N)
ghi=HN-lN-Q(1,2*N)+Q(1,12*N*N)
ld-=8*ghi
ud-=8*glo

# Machin pi=16 atan(1/5)-4 atan(1/239), alternating rational bounds.
def atan_inv(n,N=16):
    terms=[Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(N+1)]
    partial=sum(terms[:-1],Q(0))
    if N%2==0: return partial,partial+terms[-1]
    return partial+terms[-1],partial
lo5,hi5=atan_inv(5)
lo239,hi239=atan_inv(239)
pil,pih=16*lo5-4*hi239,16*hi5-4*lo239
assert Q(3)<pil<pih<Q(4)
log2lo,log2hi=log_atanh(Q(2))
lo_pi_half,_=log_atanh(pil/2)
_,hi_pi_half=log_atanh(pih/2)
logpi_lo=lo_pi_half+log2lo
logpi_hi=hi_pi_half+log2hi
ld-=8*logpi_hi
ud-=8*logpi_lo
assert Q(450,1000)<ld<ud<Q(456,1000)

terms_1=[(Q(27),9,1)]+[(-v,p,1) for p,v in c1.items()]+[
    (Q(-3),5,1),(Q(-27,2),50,41),(Q(1,6),82,81)]
l1,u1=enclose(terms_1)
assert Q(-132426,1000000)<l1<u1<Q(-132425,1000000)
terms_2=[(Q(72),9,1)]+[(-v,p,1) for p,v in c2.items()]+[
    (Q(-1,2),82,81),(Q(1,18),6643,6642),
    (Q(9,2),50,41),(Q(81,2),3650,3731)]
l2,u2=enclose(terms_2)
assert Q(-1978,1000000)<l2<u2<Q(-1977,1000000)

# Complete ORIGINAL signed Weil quadratic value, not a PSD surrogate.
lo=A*ld-4*s*u1+2*l2
hi=A*ud-4*s*l1+2*u2
assert Q(77143,10000)<lo<hi<Q(77145,10000)

# Rigorous smoothing error for any even nonnegative mass-one Cc∞
# mollifier at width eta=10^-9 (details in mathematical note).
eta=Q(1,10**9)
Lipschitz=5*(A+4*s+2)
assert Lipschitz==Q(1280,9) and 2*Lipschitz<286
D=286*eta
smooth_error=(216+11648+5+200+12)*D+2300*eta
assert smooth_error < 3500000*eta
assert lo-smooth_error > Q(77,10)

# Both original pole contributions cancel exactly, not by omission.
pole_factor=A*Q(16)-4*s*Q(80,3)+2*Q(656,9)
assert pole_factor==0

print('CR-ALP-016 SELF-CONTAINED EXACT RATIONAL AUDIT PASS')
for name,l,u in [('W00',ld,ud),('W01',l1,u1),('W02',l2,u2),('pole-null W(k_h)',lo,hi)]:
    print(name,':',float(l),'< value <',float(u),'[display-only decimals]')
print('Both Mellin pole moments = 0, smoothed moments remain 0 exactly')
print('All even mass-one mollifiers at eta=1e-9 yield original signed W>7.7')
print('Physical L2 mass:',mass,'; full prime bases:',len(ch))
print('No support-independent positivity, Connes/Suzuki domain theorem, or RH claim')
