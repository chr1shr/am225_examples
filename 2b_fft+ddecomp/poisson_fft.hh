#ifndef POISSON_FFT
#define POISSON_FFT

#include "file_output.hh"
#include <fftw3.h>

class poisson_fft {
    public:
        /** The number of gridpoints in one dimension. */
        const int n;
        /** The total number of gridpoints. */
        const int nn;
        /** The grid spacing. */
        const double h;
        /** The discretized source term in the Poisson equation. */
        double* const f;
        /** The discretized solution to the Poisson equation. */
        double* const v;
        /** The frequency domain. */
        double* const w;
        poisson_fft(int n_);
        ~poisson_fft();
        void init();
        void solve();
        void init_mms();
        double l2_error_mms();
        void print(bool solution);
        void output_solution(const char* filename);
        /** Outputs the source term in the 2D Gnuplot matrix format.
         * \param[in] filename the name of the file to write to. */
        inline void output_source(const char* filename) {
            gnuplot_output(filename,f,n,n,h,1-h,h,1-h);
        }
    private:
        /** An array holding the eigenvalues of the one-dimensional Poisson
         * matrix T_N. */
        double* const lam;
        /** The FFTW plan for converting the source term into the frequency
         * domain. */
        fftw_plan plan_fwd;
        /** The FFTW plan for converting the frequency domain back to the
         * solution. */
        fftw_plan plan_bck;
};

#endif
