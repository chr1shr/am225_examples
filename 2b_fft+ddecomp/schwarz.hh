#ifndef SCHWARZ_HH
#define SCHWARZ_HH

#include "poisson_fft.hh"

class schwarz {
    public:
        /** The size of the square subdomains. */
        const int n;
        /** The x offset of the second square subdomain. */
        const int xo;
        /** The y offset of the second square subdomain. */
        const int yo;
        /** The total x-dimension of the augmented grid. */
        const int ma;
        /** The total y-dimension of the augmented grid. */
        const int na;
        /** The total number of grid points in the augmented grid. */
        const int mna;
        /** The solution on the augmented grid. */
        double* const v;
        /** The source term on the augmented grid. */
        double* const f;
        /** The residual on the augmented grid. */
        double* const r;
        schwarz(int n_,int xo_,int yo_);
        ~schwarz();
        void init();
        void solve_additive(bool verbose=true);
        void solve_multiplicative(bool verbose=true);
        void clear_solution();
        inline void output_solution(const char* filename) {
            gnuplot_output(filename,v,ma,na,-1,1+grid1.h*xo,-1,1+grid1.h*yo);
        }
        inline void output_source(const char* filename) {
            gnuplot_output(filename,f,ma,na,-1,1+grid1.h*xo,-1,1+grid1.h*yo);
        }
        inline void output_residual(const char* filename) {
            gnuplot_output(filename,r,ma,na,-1,1+grid1.h*xo,-1,1+grid1.h*yo);
        }
    private:
        double compute_residual();
        inline double store_resid(int k);
        inline double resid(int k);
        double residual_to_square(int o,double *g);
        void copy_to_square(double *p,double *g);
        void add_from_square(double *g,double *p);
        /** The first square subdomain. */
        poisson_fft grid1;
        /** The second square subdomain. */
        poisson_fft grid2;
        /** The tolerance on a line in the linear system, taking into account
         * the machine precision, the scale of an equation in the linear system
         * (roughly set by 1/h^2) and a small multiplicative padding factor. */
        const double tol;
        /** The square tolerance for the additive method, normalized by the
         * total number of gridpoints. */
        const double ntsq_a;
        /** The square tolerance for the multiplicative method, normalized by
         * the total number of gridpoints. */
        const double ntsq_m;
};

#endif
