#include "msu_fem.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
multisetup1::multisetup1(fluid_2d &f) : m(f.m), n(f.n), mn(f.mn), x_prd(f.x_prd),
    y_prd(f.y_prd), acc(1e-20), xxsp(f.xxsp), yysp(f.yysp), z(f.smac) {}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *              respectively.
 * \param[in] z_ a pointer to the solution array. */
multisetup1::multisetup1(int m_,int n_,bool x_prd_,bool y_prd_,
        double dx,double dy,double *z_) : m(m_),
    n(n_), mn(m*n), x_prd(x_prd_), y_prd(y_prd_), acc(1e-20),
    xxsp(1./(dx*dx)), yysp(1./(dy*dy)), z(z_) {}

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
multisetup2::multisetup2(fluid_2d &f) : m(f.m_fem), n(f.n_fem), mn(m*n),
    x_prd(f.x_prd), y_prd(f.y_prd), acc(1e-20), dydx(f.dy/f.dx),
    dxdy(f.dx/f.dy), fm(4./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-1./6.*(dxdy+dydx)), z(f.sfem) {}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *              respectively.
 * \param[in] z_ a pointer to the solution array. */
multisetup2::multisetup2(int m_,int n_,bool x_prd_,bool y_prd_,
        double dx,double dy,double *z_) : m(m_),
    n(n_), mn(m*n), x_prd(x_prd_), y_prd(x_prd_), acc(1e-20), dydx(dy/dx),
    dxdy(dx/dy), fm(4./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-1./6.*(dxdy+dydx)), z(z_) {}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<multisetup1,double,double>;
template void tgmg_base<multisetup1,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<multisetup1,double,double>::clear_z();
template void tgmg_base<multisetup1,double,double>::output_res(char const*,double,double,double,double);
template class tgmg<multisetup2,double,double>;
template void tgmg_base<multisetup2,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<multisetup2,double,double>::output_res(char const*,double,double,double,double);
template void tgmg_base<multisetup2,double,double>::clear_z();
template class tgmg_base<tgmg_level<double,double>,double,double>;
