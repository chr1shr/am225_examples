#ifndef MGS_FEM_VR_HH
#define MGS_FEM_VR_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "tgmg.hh"

class fluid_2d;

struct mgs_fem_varying_rho {
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
    /** The Gauss-Seidel mode to use. */
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
    /** An array for holding the source term, if the class is initialized
     * in stand-alone mode. */
    double* const src;
    /** An array for holding the reciprocal of the density field. */
    double* const irho;
    enum cell_center {dl,dr,ul,ur};
    mgs_fem_varying_rho(fluid_2d &f);
    mgs_fem_varying_rho(int m_,int n_,bool x_prd_,bool y_prd_,double dx,double dy);
    ~mgs_fem_varying_rho();
    inline bool not_l(int i) {return x_prd||i>0;}
    inline bool not_r(int i) {return x_prd||i<m-1;}
    inline bool not_lr(int i) {return x_prd||(i>0&&i<m-1);}
    inline bool not_d(int ij) {return y_prd||ij>=m;}
    inline bool not_u(int ij) {return y_prd||ij<mn-m;}
    inline bool not_du(int ij) {return y_prd||(ij>=m&&ij<mn-m);}
    inline bool not_dl(int i,int ij) {return not_d(ij)&&not_l(i);}
    inline bool not_dr(int i,int ij) {return not_d(ij)&&not_r(i);}
    inline bool not_ul(int i,int ij) {return not_u(ij)&&not_l(i);}
    inline bool not_ur(int i,int ij) {return not_u(ij)&&not_r(i);}
    inline double cp(int i,int ij,int dir) {
        switch(dir) {
            case dl: return irho[ij];
            case dr: return x_prd&&i==m-1?irho[ij+1-m]:irho[ij+1];
            case ul: return y_prd&&ij>=mn-m?irho[ij+m-mn]:irho[ij+m];
        }
        return irho[ij+(x_prd&&i==m-1?1-m:1)+(y_prd&&ij>=mn-m?m-mn:m)];
    }
    inline double a_dl(int i,int ij) {return not_dl(i,ij)?fc*cp(i,ij,dl):0;}
    inline double a_dr(int i,int ij) {return not_dr(i,ij)?fc*cp(i,ij,dr):0;}
    inline double a_ul(int i,int ij) {return not_ul(i,ij)?fc*cp(i,ij,ul):0;}
    inline double a_ur(int i,int ij) {return not_ur(i,ij)?fc*cp(i,ij,ur):0;}
    inline double a_dc(int i,int ij) {
        const double rho_dc=(not_l(i)?cp(i,ij,dl):0) + (not_r(i)?cp(i,ij,dr):0);
        return not_d(ij)?rho_dc*hey:0;
    }
    inline double a_uc(int i,int ij) {
        const double rho_uc=(not_l(i)?cp(i,ij,ul):0) + (not_r(i)?cp(i,ij,ur):0);
        return not_u(ij)?rho_uc*hey:0;
    }
    inline double a_cl(int i,int ij) {
        const double rho_cl=(not_d(ij)?cp(i,ij,dl):0) + (not_u(ij)?cp(i,ij,ul):0);
        return not_l(i)?rho_cl*hex:0;
    }
    inline double a_cr(int i,int ij) {
        const double rho_cr=(not_d(ij)?cp(i,ij,dr):0) + (not_u(ij)?cp(i,ij,ur):0);
        return not_r(i)?rho_cr*hex:0;
    }
    inline double a_cc(int i,int ij) {
        return (!not_d(ij)?((not_l(i)?cp(i,ij,ul):0)+(not_r(i)?cp(i,ij,ur):0))
              :(!not_u(ij))?((not_l(i)?cp(i,ij,dl):0)+(not_r(i)?cp(i,ij,dr):0))
              :((not_l(i)?(cp(i,ij,ul)+cp(i,ij,dl)):0)
               +(not_r(i)?(cp(i,ij,ur)+cp(i,ij,dr)):0)))*fm;
    }
    inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
    double mul_a(int i,int ij);
    /** Sets up the linear systems on the coarser grids. */
    inline void setup() {mg.setup();}
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
    tgmg<mgs_fem_varying_rho,double,double> mg;
};

#endif
