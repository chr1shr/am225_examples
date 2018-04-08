#include "mgs_fem.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_fem::mgs_fem(fluid_2d &f) : m(f.m_fem), n(f.n_fem), mn(m*n),
    x_prd(f.x_prd), y_prd(f.y_prd), dydx(f.dy/f.dx), dxdy(f.dx/f.dy),
    fm(4./3.*(dxdy+dydx)), fm_inv(1.0/fm), fey(1./3.*(-2*dxdy+dydx)),
    hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)), hex(0.5*fex),
    fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(fm,1e4)), z(new double[mn]),
    mg(*this,f.src,z) {
    mg.setup();
    mg.clear_z();
}

/** Evaluates a single entry of the linear matrix system multiplied by the
 * current vector (held in the z array).
 * \param[in] i the horizontal index of the entry to consider.
 * \param[in] ij the index of the entry to consider.
 * \return The component of the matrix-vector product. */
double mgs_fem::mul_a(int i,int ij) {
    double *w=z+ij;

    // Since most gridpoints are interior, deal with this case
    // first
    if(i>0&&i<m-1&&ij>=m&&ij<mn-m) return fc*(w[-m-1]+w[-m+1]
            +w[m-1]+w[m+1])+fey*(w[-m]+w[m])+fex*(w[-1]+w[1]);

    // Compute constants and memory shifts for x-periodicity
    int sl=-1,sr=1;
    double lmu=1,rmu=1,afex=fex,afey=fey,ans;
    if(i==0) {
        sl+=m;
        if(!x_prd) {lmu=0;afey=hey;}
    } else if(i==m-1) {
        sr-=m;
        if(!x_prd) {rmu=0;afey=hey;}
    }

    // Assemble terms, taking into account y-periodicity
    if(ij>=m) ans=fc*(lmu*w[-m+sl]+rmu*w[-m+sr])+afey*w[-m];
    else if(y_prd) ans=fc*(lmu*w[mn-m+sl]+rmu*w[mn-m+sr])+afey*w[mn-m];
    else {ans=0;afex=hex;}
    if(ij<mn-m) ans+=fc*(lmu*w[m+sl]+rmu*w[m+sr])+afey*w[m];
    else if(y_prd) ans+=fc*(lmu*w[m-mn+sl]+rmu*w[m-mn+sr])+afey*w[m-mn];
    else afex=hex;
    return ans+afex*(lmu*w[sl]+rmu*w[sr]);
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_fem,double,double>;
template void tgmg_base<mgs_fem,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_fem,double,double>::output_res(char const*,double,double,double,double);
template void tgmg_base<mgs_fem,double,double>::clear_z();
template class tgmg_base<tgmg_level<double,double>,double,double>;
