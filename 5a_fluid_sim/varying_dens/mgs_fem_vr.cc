#include "mgs_fem_vr.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_fem_varying_rho::mgs_fem_varying_rho(fluid_2d &f) : m(f.m_fem), n(f.n_fem),
    mn(m*n), x_prd(f.x_prd), y_prd(f.y_prd), dydx(f.dy/f.dx), dxdy(f.dx/f.dy),
    fm(1./3.*(dxdy+dydx)), fm_inv(1.0/fm), fey(1./3.*(-2*dxdy+dydx)),
    hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)), hex(0.5*fex),
    fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(fm,1e4)), z(new double[mn]),
    src(NULL), irho(new double[mn]), mg(*this,f.src,z) {
    mg.clear_z();
}

mgs_fem_varying_rho::mgs_fem_varying_rho(int m_,int n_,bool x_prd_,bool y_prd_,double dx,double dy) :
    m(m_), n(n_), mn(m*n), x_prd(x_prd_), y_prd(y_prd_), dydx(dy/dx),
    dxdy(dx/dy), fm(1./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(fm,1e4)),
    z(new double[mn]), src(new double[mn]), irho(new double[mn]), mg(*this,src,z) {
    mg.clear_z();
}

mgs_fem_varying_rho::~mgs_fem_varying_rho() {
    delete [] irho;
    if(src!=NULL) delete [] src;
    delete [] z;
}

/** Evaluates a single entry of the linear matrix system multiplied by the
 * current vector (held in the z array).
 * \param[in] i the horizontal index of the entry to consider.
 * \param[in] ij the index of the entry to consider.
 * \return The component of the matrix-vector product. */
double mgs_fem_varying_rho::mul_a(int i,int ij) {
    double *w=z+ij;
    const double cul=cp(i,ij,ul),cdl=cp(i,ij,dl),
                 cur=cp(i,ij,ur),cdr=cp(i,ij,dr);

    // Since most gridpoints are interior, deal with this case
    // first
    if(i>0&&i<m-1&&ij >= m && ij < mn-m)
    return cul*(w[-1]*hex+w[+m-1]*fc+w[+m]*hey)
          +cdl*(w[-1]*hex+w[-m-1]*fc+w[-m]*hey)
          +cur*(w[+1]*hex+w[+m+1]*fc+w[+m]*hey)
          +cdr*(w[+1]*hex+w[-m+1]*fc+w[-m]*hey);

    // Compute constants and memory shifts for x-periodicity
    int sl=-1,sr=1;
    double lmu=1,rmu=1,ans,hex_l=hex*(cdl+cul),hex_r=hex*(cdr+cur),
         hey_d=hey*(cdl+cdr),hey_u=hey*(cul+cur);
    if(i==0) {
        sl+=m;
        if(!x_prd) {lmu=0;hey_d=hey*cdr;hey_u=hey*cur;}
    } else if(i==m-1) {
        sr-=m;
        if(!x_prd) {rmu=0;hey_d=hey*cdl;hey_u=hey*cul;}
    }

    // Assemble terms,taking into account y-periodicity
    if(ij>=m) ans=fc*(lmu*w[-m+sl]*cdl+rmu*w[-m+sr]*cdr)+hey_d*w[-m];
    else if(y_prd) ans=fc*(lmu*w[mn-m+sl]*cdl+rmu*w[mn-m+sr]*cdr)+hey_d*w[mn-m];
    else ans=0,hex_l=hex*cul,hex_r=hex*cur;
    if(ij<mn-m) ans+=fc*(lmu*w[m+sl]*cul+rmu*w[m+sr]*cur)+hey_u*w[m];
    else if(y_prd) ans+=fc*(lmu*w[m-mn+sl]*cul+rmu*w[m-mn+sr]*cur)+hey_u*w[m-mn];
    else hex_l=hex*cdl,hex_r=hex*cdr;
    return ans+hex_l*lmu*w[sl]+hex_r*rmu*w[sr];
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_fem_varying_rho,double,double>;
template void tgmg_base<mgs_fem_varying_rho,double,double>
    ::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_fem_varying_rho,double,double>
    ::clear_z();
template void tgmg_base<mgs_fem_varying_rho,double,double>
    ::output_res(char const*,double,double,double,double);
template double tgmg_base<mgs_fem_varying_rho,double,double>
    ::mds();
template class tgmg_base<tgmg_level<double,double>,double,double>;
