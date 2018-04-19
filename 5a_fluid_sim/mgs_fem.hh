#ifndef MGS_FEM_HH
#define MGS_FEM_HH

class fluid_2d;

#include "tgmg.hh"

struct mgs_fem {
    /** The number of gridpoints in the x direction. */
    const int m;
    /** The number of gridpoints in the y direction. */
    const int n;
    /** Total number of gridpoints. */
    const int mn;
    /** Periodicity in the x direction. */
    const bool x_prd;
    /** Periodicity in the y direction. */
    const bool y_prd;
    /** The mode to use for the Gauss-Seidel smoothing. (0=default) */
    static const char gs_mode=0;
    /** The ratio of vertical and horizontal grid spacings. */
    const double dydx;
    /** The ratio of horizontal and vertical grid spacings. */
    const double dxdy;
    /** The central term in the finite element stencil. */
    const double fm;
    /** The reciprocal of the central term in the finite element stencil,
     * needed during the Gauss--Seidel smoothing steps. */
    const double fm_inv;
    /** The vertical term in the finite element stencil. */
    const double fey;
    /** Half the vertical term in the finite element stencil, required
     * at boundary gridpoints. */
    const double hey;
    /** The horizontal term in the finite element stencil. */
    const double fex;
    /** Half the horizontal term in the finite element stencil, required
     * at boundary gridpoints. */
    const double hex;
    /** The diagonal corner term in the finite element stencil. */
    const double fc;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** The array holding the solution of the linear system (i.e. the fluid
     * pressure). */
    double* const z;
    mgs_fem(fluid_2d &f);
    ~mgs_fem() {
        delete [] z;
    }
    inline bool not_l(int i) {return x_prd||i>0;}
    inline bool not_r(int i) {return x_prd||i<m-1;}
    inline bool not_lr(int i) {return x_prd||(i>0&&i<m-1);}
    inline bool not_d(int ij) {return y_prd||ij>=m;}
    inline bool not_u(int ij) {return y_prd||ij<mn-m;}
    inline bool not_du(int ij) {return y_prd||(ij>=m&&ij<mn-m);}
    inline double a_dl(int i,int ij) {return not_d(ij)&&not_l(i)?fc:0;}
    inline double a_dr(int i,int ij) {return not_d(ij)&&not_r(i)?fc:0;}
    inline double a_ul(int i,int ij) {return not_u(ij)&&not_l(i)?fc:0;}
    inline double a_ur(int i,int ij) {return not_u(ij)&&not_r(i)?fc:0;}
    inline double a_dc(int i,int ij) {return not_d(ij)?(not_lr(i)?fey:hey):0;}
    inline double a_uc(int i,int ij) {return not_u(ij)?(not_lr(i)?fey:hey):0;}
    inline double a_cl(int i,int ij) {return not_l(i)?(not_du(ij)?fex:hex):0;}
    inline double a_cr(int i,int ij) {return not_r(i)?(not_du(ij)?fex:hex):0;}
    inline double a_cc(int i,int ij) {
        return (not_lr(i)?1:0.5)*(not_du(ij)?1:0.5)*fm;
    }
    inline double inv_cc(int i,int ij,double v) {
        return (not_lr(i)?1:2)*(not_du(ij)?1:2)*fm_inv*v;
    }
    double mul_a(int i,int ij);
    /** Solves the linear system using multigrid V-cycles. */
    inline void solve_v_cycle() {
        if(!mg.solve_v_cycle(tp)) {
            fputs("V-cycle failed to converge in FEM problem\n",stderr);
            exit(1);
        }
    }
    /** A helper class for the multigrid library that holds information for
     * predicting the number of V-cycles that are required. */
    tgmg_predict tp;
    /** The multigrid solver. */
    tgmg<mgs_fem,double,double> mg;
};

#endif
