#ifndef HO_TRANSPORT_HH
#define HO_TRANSPORT_HH

#include <cstdio>
#include <cmath>

class traffic {
    public:
        /** The grid size. */
        const int m;
        /** The grid spacing. */
        const double dx;
        /** The timestep. */
        double dt;
        /** The primary grid for storing the solution. */
        double *a;
        /** The secondary grid for storing the solution. */
        double *b;
        /** An array for the fluxes. */
        double* const F;
        traffic(int m_,double A_);
        ~traffic();
	void godunov(double dt);
        void init_step_function();
        void init_exp_sine();
        void solve(const char* filename,int snaps,double duration,double safe_fac);
        double integral();
    private:
	inline double flux(double x) {
		return x*(1-x);
	}
        void print_line(FILE *fp,double x,double *zp,int snaps);
};

#endif
