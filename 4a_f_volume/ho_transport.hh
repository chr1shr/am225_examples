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
        ho_transport(int m_,double A_);
        ~ho_transport();
        void lax_wendroff(double dt);
        void beam_warming(double dt);
        void init_step_function();
        void init_exp_sine();
        void init_nu_array(int type,double safe_fac);
        void solve(const char* filename,int snaps,double duration,double safe_fac,int type);
        double integral();
    private:
        void print_line(FILE *fp,double x,double *zp,int snaps);
};

#endif
