// Global second projected signed-prime moment; includes the frozen 058 evaluator.
// Build against the extracted, SHA-verified original source; no replacement coefficients.
#define main cr_alp_058_original_main
#include "validated.cpp"
#undef main
struct Move59{int j,s,p;};
struct Panel59{I m,h;vector<Move59> active;};
I sq59(I x){Q z=x.absmax(),l=x.l<=0&&x.h>=0?0:min(fabsq(x.l),fabsq(x.h));return I(l==0?0:down(l*l),up(z*z));}
int main(int argc,char**argv){
 assert(argc==5&&fegetround()==FE_TONEAREST);
 mpfr_init2(mm,192);mpfr_init2(nn,192);mpfr_init2(oo,192);
 init(argv[1],argv[2],argv[3]);
 // Same first Q85 projection as CR-ALP-058. Reintegrate, not a fitted finite projection.
 vector<I> proj(nl);I vb;
 for(auto &p:mesh[0])for(int j=0;j<96;j++){
  I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];vector<I>b(nl);I x=eval(t,p.p,&b);
  for(int i=0;i<nl;i++)proj[i]=proj[i]+w*x*b[i];vb=vb+w*x*bval(t);
 }
 vb=widen(vb,1e-12Q);I st=(ex(a)-ex(-a))/I(2)+I(sgn)*a;vector<I> c(nl);
 for(int i=0;i<nl;i++){
  proj[i]=widen(proj[i],1e-12Q);I v=nu*I(ids[i]);
  c[i]=I(2)*(bp+bm)*v/(rt*(v*v+I(1)/I(4)));st=st-sq59(c[i]);vb=vb-c[i]*proj[i];
 }
 assert(st.l>0.003Q);I lam=vb/st;
 for(auto &p:pieces){for(int i=0;i<nl;i++)p.s[i]=p.s[i]+(c[i]*lam-proj[i])/rt;p.ep=p.ep-lam*bp;p.en=p.en-lam*bm;}
 midpointize();
 I G;
 for(auto &p:mesh[0])for(int j=0;j<96;j++){I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];G=G+w*sq59(eval(t,p.p));}
 G=widen(G,1e-12Q);assert(G.l>0&&G.h<0.16Q);out("proxy_G",G);
 ifstream in(argv[4]);int n;in>>n;assert(in&&n>0);vector<Panel59> pan;
 for(int k=0;k<n;k++){Panel59 p;p.m=rd(in);p.h=rd(in);int m;in>>m;assert(m>=0&&m<=14);
  for(int j=0;j<m;j++){Move59 z;in>>z.j>>z.s>>z.p;assert(in&&z.j>=0&&z.j<7&&(z.s==1||z.s==-1)&&z.p>=0&&z.p<(int)pieces.size());p.active.push_back(z);}pan.push_back(p);}
 string trailing;assert(!(in>>trailing));
 I V2, bpv;vector<I> vproj(nl);
 for(auto &p:pan)for(int j=0;j<96;j++){
  I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j],v;
  for(auto z:p.active)v=v-wt[z.j]*eval(t+I(z.s)*sh[z.j],z.p);
  V2=V2+w*sq59(v);bpv=bpv+w*v*bval(t);
  // Direct MPFR interval sine; no new midpoint error argument.
  for(int i=0;i<nl;i++){auto sc=sincosi(nu*I(ids[i])*t);vproj[i]=vproj[i]+w*v*sc.first/rt;}
 }
 V2=widen(V2,1e-12Q);bpv=widen(bpv,1e-12Q);
 I sine2;for(int i=0;i<nl;i++){vproj[i]=widen(vproj[i],1e-12Q);sine2=sine2+sq59(vproj[i]);bpv=bpv-c[i]*vproj[i];}
 I bridge2=sq59(bpv)/st;
 I N=V2-sine2-bridge2;assert(N.h>=0);N.l=max(N.l,(Q)0);
 assert(V2.h<2&&sine2.h<2&&bridge2.h<2&&N.h<2);
 out("V2",V2);out("sine2",sine2);out("bridge2",bridge2);out("N_proxy",N);
 mpfr_clear(mm);mpfr_clear(nn);mpfr_clear(oo);
}
