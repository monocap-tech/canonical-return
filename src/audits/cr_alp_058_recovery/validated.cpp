// Scalar physical overlap enclosure, with outward binary128 intervals and a
// separately bounded midpoint evaluator. See VALIDATION.md for the forward error proof.
// MPFR C ABI declarations (libmpfr 4.x, GNU/Linux x86-64); no development header installed.
#include <quadmath.h>
#include <bits/stdc++.h>
#include <cfenv>
#include <gmp.h>
using namespace std;using Q=__float128;
extern "C" {
 typedef long mpfr_prec_t; typedef long mpfr_exp_t; typedef int mpfr_sign_t;
 typedef struct {mpfr_prec_t _mpfr_prec;mpfr_sign_t _mpfr_sign;mpfr_exp_t _mpfr_exp;mp_limb_t *_mpfr_d;} __mpfr_struct;
 typedef __mpfr_struct mpfr_t[1];
 void mpfr_init2(mpfr_t,mpfr_prec_t);void mpfr_clear(mpfr_t);
 int mpfr_set_str(mpfr_t,const char*,int,int);int mpfr_set_float128(mpfr_t,Q,int);Q mpfr_get_float128(const mpfr_t,int);
 int mpfr_const_pi(mpfr_t,int);int mpfr_sin(mpfr_t,const mpfr_t,int);int mpfr_cos(mpfr_t,const mpfr_t,int);int mpfr_exp(mpfr_t,const mpfr_t,int);int mpfr_sqrt(mpfr_t,const mpfr_t,int);
 int mpfr_sin_cos(mpfr_t,mpfr_t,const mpfr_t,int);
}
constexpr int RN=0,RU=2,RD=3;
mpfr_t mm,nn,oo;
Q mpstr(const string&s,int r){assert(mpfr_set_str(mm,s.c_str(),10,r)==0);return mpfr_get_float128(mm,r);}
Q down(Q x){return nextafterq(x,-HUGE_VALQ);} Q up(Q x){return nextafterq(x,HUGE_VALQ);}
struct I {Q l,h; I():l(0),h(0){} I(int n):l(n),h(n){} explicit I(Q n):l(n),h(n){} I(Q a,Q b):l(a),h(b){assert(a<=b);} Q mid()const{return (l+h)/2;} Q rad()const{Q m=mid();return max(up(m-l),up(h-m));} Q absmax()const{return max(fabsq(l),fabsq(h));}};
I operator+(I a,I b){return I(down(a.l+b.l),up(a.h+b.h));} I operator-(I a,I b){return I(down(a.l-b.h),up(a.h-b.l));} I operator-(I a){return I(-a.h,-a.l);} I operator*(I a,I b){Q p[]={a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h};return I(down(*min_element(p,p+4)),up(*max_element(p,p+4)));} I operator/(I a,I b){assert(b.l*b.h>0);return a*I(down(1/b.h),up(1/b.l));}
I widen(I x,Q r){return I(down(x.l-r),up(x.h+r));} I rd(istream&f){string l,h;f>>l>>h;assert(f);return I(mpstr(l,RD),mpstr(h,RU));}
Q fun(Q x,int which,int r=RN){mpfr_set_float128(mm,x,RN);if(which==0)mpfr_exp(nn,mm,r);else if(which==1)mpfr_sin(nn,mm,r);else if(which==2)mpfr_cos(nn,mm,r);else mpfr_sqrt(nn,mm,r);return mpfr_get_float128(nn,r);}
I ex(I x){return I(fun(x.l,0,RD),fun(x.h,0,RU));} I sqrti(I x){assert(x.l>0);return I(fun(x.l,3,RD),fun(x.h,3,RU));}
pair<I,I> sincosi(I x){Q m=x.mid(),r=x.rad();return {widen(I(fun(m,1,RD),fun(m,1,RU)),r),widen(I(fun(m,2,RD),fun(m,2,RU)),r)};}
void sc(Q x,Q&s,Q&c){mpfr_set_float128(mm,x,RN);mpfr_sin_cos(nn,oo,mm,RN);s=mpfr_get_float128(nn,RN);c=mpfr_get_float128(oo,RN);}
Q expn(Q x){return fun(x,0);}
string ss(Q x){char b[100];quadmath_snprintf(b,sizeof(b),"%.23Qe",x);return b;} void out(string n,I v){ // rounded display only, outward printed decimal bounds at deliberately loose 1e-15 grid
 Q t=1e15Q;long long l=(long long)floorq(down(v.l*t)),h=(long long)ceilq(up(v.h*t));cout<<n<<" "<<l<<" "<<h<<" / 1e15\n"<<flush;}
struct G{I c;vector<I>A;};struct PG{Q c;vector<Q>A;};struct Piece{vector<I>s,c;I ep,en;};struct PP{vector<Q>s,c;Q ep,en,r;};struct Panel{I m,h;int p,r;};
I a,d,ga,nu,rt,bp,bm,f,db,eb;int k,sgn,nl;vector<int> ids;vector<I>src,tr,lo,sh,wt,gn,gw;vector<G> gs;vector<PG> pg;vector<Piece>pieces;vector<PP>pp;vector<vector<Panel>> mesh;Q NUF,RT,BP,BM,F,EB,D=2.2Q;vector<Q>plo;Q baseradius;
I bval(I t){return bp*ex(t/I(2))+bm*ex(-t/I(2));}
Q group(Q t,const PG&g){Q x=g.c*t;if(x>1000)return 0;assert(x>=0);Q v=expn(-x),sum=0;for(int j=0;j<(int)g.A.size();j++){if(j)v*=x/j;sum+=v*g.A[j];}return sum;}
Q boundary(Q t){Q v=expn(-t/2),w=expn(-(D-t)/2),u=expn(-2*t),z=expn(-2*(D-t));Q val=0;for(Q A:plo){val+=A*(v+sgn*w);v*=u;w*=z;}for(auto&g:pg)val+=group(t,g)+sgn*group(D-t,g);Q et=expn(t/2),em=1/et;return val-F*(BP*(et*(1+D-t)-em)+BM*(em*(t+1)-et*EB));}
void midpointize(){pp.clear();for(auto &p:pieces){PP v;v.ep=p.ep.mid();v.en=p.en.mid();I r=I(4)*(I(p.ep.rad())+I(p.en.rad())),mag=I(0);for(int i=0;i<(int)ids.size();i++){v.s.push_back(p.s[i].mid());v.c.push_back(p.c[i].mid());r=r+I(p.s[i].rad())+I(p.c[i].rad());mag=mag+I(p.s[i].absmax())+I(p.c[i].absmax());}assert(mag.h<1e17Q);assert(p.ep.absmax()<1e17Q&&p.en.absmax()<1e17Q);v.r=(r+I(baseradius)+I(1e-8Q)).h;assert(v.r<0.01Q);pp.push_back(v);}Q mx=0;for(auto &p:pp)mx=max(mx,p.r);cerr<<"max point radius "<<ss(mx)<<"\n";}
I eval(I ti,int p,vector<I>*basis=nullptr){assert(ti.rad()<1e-31Q&&ti.l>=0&&ti.h<=up(D));Q t=ti.mid();auto&v=pp[p];Q sum=boundary(t)+v.ep*expn(t/2)+v.en*expn(-t/2);Q sa,ca,c=0,s=0;sc(2*NUF*t,sa,ca);
 for(int j=0;j<(int)ids.size();j++){if(j%16==0)sc(NUF*ids[j]*t,s,c);else{Q next=s*ca+c*sa;c=c*ca-s*sa;s=next;}sum+=v.s[j]*s+v.c[j]*c;if(basis&&j<nl)(*basis)[j]=widen(I(s/RT),1e-26Q);}
 return widen(I(sum),v.r);}
void init(const char*input,const char*meshfile,const char*nodes){
 a=I(11)/I(10);d=I(2)*a;ga=I(72)/I(25);mpfr_const_pi(mm,RD);Q pil=mpfr_get_float128(mm,RD);mpfr_const_pi(mm,RU);Q pih=mpfr_get_float128(mm,RU);nu=I(pil,pih)/d;rt=sqrti(a);NUF=nu.mid();RT=rt.mid();assert(nu.rad()<1e-32Q&&rt.rad()<1e-32Q&&d.rad()<1e-32Q);
 ifstream in(input);int n;in>>k;sgn=k==1?1:-1;f=rd(in);db=rd(in);in>>n;nl=0;for(int j=0;j<n;j++){int i;in>>i;ids.push_back(i);if(i<=85)nl++;src.push_back(rd(in));tr.push_back(rd(in));}in>>n;for(int j=0;j<n;j++){rd(in);lo.push_back(rd(in));}in>>n;for(int j=0;j<n;j++){G g;g.c=rd(in);int nr;in>>nr;for(int r=0;r<nr;r++)g.A.push_back(rd(in));gs.push_back(g);}
 bp=ex(-a/I(2))/I(2);bm=I(sgn)*ex(a/I(2))/I(2);eb=ex(-d);BP=bp.mid();BM=bm.mid();F=f.mid();EB=eb.mid();assert(bp.rad()<1e-32Q&&bm.rad()<1e-32Q&&eb.rad()<1e-32Q&&f.absmax()<100);I r=I(0),lm=I(0),gm=I(0);for(auto x:lo){plo.push_back(x.mid());r=r+I(2)*I(x.rad());lm=lm+I(x.absmax());}for(auto &g:gs){PG h;h.c=g.c.mid();assert(g.c.rad()<1e-20Q&&g.c.h<17179869184.0Q);for(auto x:g.A){h.A.push_back(x.mid());r=r+I(2)*I(x.rad());gm=gm+I(x.absmax());}pg.push_back(h);}assert(lm.h<1e14Q&&gm.h<1e3Q);r=r+I(100)*(I(f.rad())+I(bp.rad())+I(bm.rad())+I(eb.rad()));baseradius=r.h;
 ifstream mf(meshfile);mf>>n;vector<vector<pair<int,int>>>active(n);for(int m=0;m<n;m++){int cnt;mf>>cnt;for(int j=0;j<cnt;j++){int i,s;mf>>i>>s;active[m].push_back({i,s});}}for(int j=0;j<7;j++){sh.push_back(rd(mf));wt.push_back(rd(mf));}for(int j=0;j<8;j++){mf>>n;vector<Panel> v;for(int m=0;m<n;m++){Panel t;t.m=rd(mf);t.h=rd(mf);mf>>t.p>>t.r;v.push_back(t);}mesh.push_back(v);}assert(mf);
 vector<vector<I>>cs(7),sn(7);for(int j=0;j<7;j++)for(int i:ids){auto [s,c]=sincosi(nu*I(i)*sh[j]);cs[j].push_back(c);sn[j].push_back(s);}
 for(auto &set:active){Piece p;for(int j=0;j<(int)ids.size();j++){p.s.push_back(-ga*tr[j]/rt);p.c.push_back(I(0));}p.ep=-ga*db*bp;p.en=-ga*db*bm;for(auto [j,s]:set){I t=I(s)*sh[j];for(int i=0;i<(int)ids.size();i++){I v=-wt[j]*(src[i]-tr[i])/rt;p.s[i]=p.s[i]+v*cs[j][i];p.c[i]=p.c[i]+I(s)*v*sn[j][i];}p.ep=p.ep-wt[j]*(f-db)*bp*ex(t/I(2));p.en=p.en-wt[j]*(f-db)*bm*ex(-t/I(2));}pieces.push_back(p);}midpointize();
 ifstream nf(nodes);nf>>n;assert(n==96);for(int j=0;j<n;j++){gn.push_back(rd(nf));gw.push_back(rd(nf));assert(gw.back().l>0);}
}
int main(int argc,char**argv){assert(argc==4);assert(fegetround()==FE_TONEAREST);mpfr_init2(mm,192);mpfr_init2(nn,192);mpfr_init2(oo,192);init(argv[1],argv[2],argv[3]);vector<I> proj(nl);I vb;
 for(auto &p:mesh[0])for(int j=0;j<96;j++){I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];vector<I>b(nl);I u=eval(t,p.p,&b);for(int i=0;i<nl;i++)proj[i]=proj[i]+w*u*b[i];vb=vb+w*u*bval(t);}
 vb=widen(vb,1e-12Q);I bm2=(ex(a)-ex(-a))/I(2)+I(sgn)*a,st=bm2;vector<I>c(nl);for(int i=0;i<nl;i++){proj[i]=widen(proj[i],1e-12Q);I v=nu*I(ids[i]);c[i]=I(2)*(bp+bm)*v/(rt*(v*v+I(1)/I(4)));st=st-c[i]*c[i];vb=vb-c[i]*proj[i];}assert(st.l>0.003Q);I lam=vb/st;cerr<<"projection lambda radius "<<ss(lam.rad())<<"\n";for(auto &p:pieces){for(int i=0;i<nl;i++)p.s[i]=p.s[i]+(c[i]*lam-proj[i])/rt;p.ep=p.ep-lam*bp;p.en=p.en-lam*bm;}midpointize();
 I G;for(auto&p:mesh[0])for(int j=0;j<96;j++){I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j],u=eval(t,p.p);G=G+w*u*u;}G=widen(G,1e-12Q);assert(G.l>0&&G.h<1);out("G",G);
 I M;for(int z=0;z<7;z++){I C;for(auto&p:mesh[z+1])for(int j=0;j<96;j++){I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];C=C+w*eval(t,p.p)*eval(t+sh[z],p.r);}C=widen(C,1e-12Q);M=M-I(2)*wt[z]*C;out("corr"+to_string(z),C);}out("M",M);I ell=I(51)/I(200),u=I(1101)/I(200);out("penalty",(ga*G-M)/(ell*u));mpfr_clear(mm);mpfr_clear(nn);mpfr_clear(oo);
}
