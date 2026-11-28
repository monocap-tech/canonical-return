// Global fourth projected signed-prime moment; first reproduce second moment; includes the frozen 058 evaluator.
// Build against the extracted, SHA-verified original source; no replacement coefficients.
#define main cr_alp_058_original_main
#include "validated.cpp"
#undef main
struct Move59{int j,s,p;};
struct Panel59{I m,h;vector<Move59> active;};
I sq59(I x){Q z=x.absmax(),l=x.l<=0&&x.h>=0?0:min(fabsq(x.l),fabsq(x.h));return I(l==0?0:down(l*l),up(z*z));}
int main(int argc,char**argv){
 assert(argc==6&&fegetround()==FE_TONEAREST);
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
 // Use exact orthonormal projection coordinates and freeze their binary128 midpoints.
 // This avoids propagating correlated projection-coefficient uncertainty pointwise.
 I sqrtst=sqrti(st),bridge=bpv/sqrtst;
 vector<Q> midcoef; I err2=sq59(I(bridge.rad()));
 I coeffsum;
 for(int i=0;i<nl;i++){
   midcoef.push_back(vproj[i].mid());err2=err2+sq59(I(vproj[i].rad()));
   coeffsum=coeffsum+I(fabsq(midcoef.back()));
 }
 Q midbridge=bridge.mid();
 assert(coeffsum.h<10&&fabsq(midbridge)<10&&st.l>0.003Q);
 // c_i are exact sine pairings enclosed by c[i]; sum magnitude bounded here.
 I csum;for(auto x:c)csum=csum+I(x.absmax());assert(csum.h<100);
 Q eta=sqrti(I(max(err2.h,1e-60Q))).h;
 assert(eta<0.1Q);out("projection_error",I(0,eta));
 // v0=A r_proxy-C, C the frozen finite correction. Thus A v0=A^2 r_proxy-A C.
 auto correction=[&](I t){
   I z,bperp=bval(t);
   for(int i=0;i<nl;i++){
     I e=sincosi(nu*I(ids[i])*t).first/rt;
     z=z+I(midcoef[i])*e;bperp=bperp-c[i]*e;
   }
   return z+I(midbridge)*bperp/sqrtst;
 };
 struct Group63{I z,weight;int p;};
 ifstream fourth(argv[5]);int n4;fourth>>n4;assert(fourth&&n4>0);
 I W2,wb;vector<I> wp(nl);
 for(int k4=0;k4<n4;k4++){
   I m=rd(fourth),h=rd(fourth);int no,ng;fourth>>no;assert(fourth&&no>=0&&no<=14);
   vector<pair<int,int>> outer;
   for(int j=0;j<no;j++){int i,s;fourth>>i>>s;assert(fourth&&i>=0&&i<7&&(s==1||s==-1));outer.push_back({i,s});}
   fourth>>ng;assert(fourth&&ng>=0&&ng<=196);vector<Group63> groups;
   for(int j=0;j<ng;j++){
     Group63 gr;gr.z=rd(fourth);int np;fourth>>gr.p>>np;assert(fourth&&gr.p>=0&&gr.p<(int)pieces.size()&&np>0&&np<=196);
     for(int k=0;k<np;k++){
       int i,s,j,t;fourth>>i>>s>>j>>t;assert(fourth&&i>=0&&i<7&&j>=0&&j<7&&(s==1||s==-1)&&(t==1||t==-1));
       I z=I(s)*sh[i]+I(t)*sh[j];assert(z.l<=gr.z.h&&gr.z.l<=z.h);
       gr.weight=gr.weight+wt[i]*wt[j];
     }
     groups.push_back(gr);
   }
   for(int j=0;j<96;j++){
     I t=m+h*gn[j],w=I(2)*h*gw[j],x;
     for(auto gr:groups)x=x+gr.weight*eval(t+gr.z,gr.p);
     for(auto [i,s]:outer)x=x+wt[i]*correction(t+I(s)*sh[i]);
     W2=W2+w*sq59(x);wb=wb+w*x*bval(t);
     for(int i=0;i<nl;i++)wp[i]=wp[i]+w*x*sincosi(nu*I(ids[i])*t).first/rt;
   }
   if((k4+1)%100==0)cerr<<"fourth panels "<<k4+1<<"/"<<n4<<"\n";
 }
 assert(!(fourth>>trailing));W2=widen(W2,1e-12Q);wb=widen(wb,1e-12Q);
 I sine4;for(int i=0;i<nl;i++){wp[i]=widen(wp[i],1e-12Q);sine4=sine4+sq59(wp[i]);wb=wb-c[i]*wp[i];}
 I bridge4=sq59(wb)/st,S=W2-sine4-bridge4;
 assert(S.h>=0&&W2.h<10&&sine4.h<10&&bridge4.h<10);S.l=max(S.l,(Q)0);
 out("W2",W2);out("sine4",sine4);out("bridge4",bridge4);out("S_proxy",S);
 mpfr_clear(mm);mpfr_clear(nn);mpfr_clear(oo);
}
