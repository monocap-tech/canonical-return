// CR-ALP-074 parity-general physical norm only.
// Include the exact authenticated 058 evaluator, unchanged. The original
// validation argument applies only if the actual input and runtime budgets pass.
#define main cr_alp_058_original_main
#include "validated.cpp"
#undef main
I sq074(I v){
    Q h=v.absmax(),l=(v.l<=0&&v.h>=0)?0:min(fabsq(v.l),fabsq(v.h));
    return I(l==0?0:down(l*l),up(h*h));
}
int main(int argc,char**argv){
    assert(argc==4&&fegetround()==FE_TONEAREST);
    mpfr_init2(mm,192);mpfr_init2(nn,192);mpfr_init2(oo,192);
    init(argv[1],argv[2],argv[3]);
    assert(k==1||k==2);
    // This is the same full Q85 projection as the original 058 scalar audit.
    vector<I> proj(nl); I vb;
    for(auto &p:mesh[0])for(int j=0;j<96;j++){
        I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];
        vector<I> b(nl);I x=eval(t,p.p,&b);
        for(int i=0;i<nl;i++)proj[i]=proj[i]+w*x*b[i];
        vb=vb+w*x*bval(t);
    }
    vb=widen(vb,1e-12Q);
    I st=(ex(a)-ex(-a))/I(2)+I(sgn)*a; vector<I> c(nl);
    for(int i=0;i<nl;i++){
        proj[i]=widen(proj[i],1e-12Q);
        I v=nu*I(ids[i]);
        c[i]=I(2)*(bp+bm)*v/(rt*(v*v+I(1)/I(4)));
        st=st-sq074(c[i]); vb=vb-c[i]*proj[i];
    }
    assert(st.l>0.003Q);
    I lam=vb/st;
    for(auto &p:pieces){
        for(int i=0;i<nl;i++)p.s[i]=p.s[i]+(c[i]*lam-proj[i])/rt;
        p.ep=p.ep-lam*bp;p.en=p.en-lam*bm;
    }
    midpointize();
    I G;
    for(auto &p:mesh[0])for(int j=0;j<96;j++){
        I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];
        G=G+w*sq074(eval(t,p.p));
    }
    G=widen(G,1e-12Q);
    // A nonnegative square does not require a strictly positive lower endpoint.
    // Reject implausible/oversized output before decimal conversion.
    assert(G.h>=0&&G.h<100);
    out("G",I(max(Q(0),G.l),G.h));
    mpfr_clear(mm);mpfr_clear(nn);mpfr_clear(oo);
}
