#ifndef DIFFUSE_HH
#define DIFFUSE_HH

#include <cstdio>
#include <cmath>

class diffuse {
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
        /** The table of integration coefficients, which incorporate the
         * timestep, the grid spacing, and the spatially dependent diffusivity.
         */
        double *nu;
        diffuse(int m_);
        ~diffuse();
        template<int type>
        void step_forward();
        void init();
        void init_nu_array(int type,double safe_fac);
        void solve(const char* filename,int snaps,int iters);
    private:
        void print_line(FILE *fp,double x,double *zp,int snaps);
        /** Calculates the diffusion constant at a position.
         * \param[in] x the position to consider.
         * \return The diffusion constant. */
        inline double f_beta(double x) {
            return 0.12+0.08*sin(2*M_PI*x);
        }
        /** The maximum of the beta function. */
        inline double beta_max() {return 0.2;}
};

#endif
