#ifndef MSU_FEM_HH
#define MSU_FEM_HH

class fluid_2d;

struct msu_fem {
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
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    const double dydx;
    const double dxdy;
    const double fm;
    const double fm_inv;
    const double fey;
    const double hey;
    const double fex;
    const double hex;
    const double fc;
    double* const z;
    multisetup2(fluid_2d &f);
    multisetup2(int m_,int n_,bool x_prd_,bool y_prd_,double dx,double dy,double *z_);
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
    inline double mul_a(int i,int ij) {
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
};

#endif
