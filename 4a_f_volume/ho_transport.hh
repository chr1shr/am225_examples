#ifndef HO_TRANSPORT_HH
#define HO_TRANSPORT_HH

#include <cstdio>
#include <cmath>

class ho_transport {
    public:
        /** The grid size. */
        const int m;
        /** The grid spacing. */
        const double dx;
        /** The advection velocity. */
        const double A;
        /** The timestep. */
        double dt;
        /** The primary grid for storing the solution. */
        double *a;
        /** The secondary grid for storing the solution. */
        double *b;
        /** An array for the slopes, sigma, multiplied by the grid spacing dx. */
        double *sdx;
        ho_transport(int m_,double A_);
        ~ho_transport();
        void godunov(double dt);
        void lax_wendroff(double dt);
        void beam_warming(double dt);
        void slope_limiter(double dt,int type);
        template<int type>
        void sl_setup();
        void step_eno2(double dt);
        void init_step_function();
        void init_exp_sine();
        void solve(const char* filename,int snaps,double duration,double safe_fac,int type);
        double integral();
    private:
        inline double eno2(double p0,double p1,double p2,double p3);
        inline double minmod(double a,double b) {
            return a*b>0?(fabs(a)<fabs(b)?a:b):0;
        }
        inline double maxmod(double a,double b) {
            return a*b>0?(fabs(a)<fabs(b)?b:a):0;
        }
        void print_line(FILE *fp,double x,double *zp,int snaps);
};

#endif
