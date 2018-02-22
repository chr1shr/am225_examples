#ifndef POISSON_FFT
#define POISSON_FFT

#include <fftw3.h>

class poisson_fft {
    public:
        const int n;
        const int nn;
        const double h;
        double* const f;
        double* const v;
        double* const w;
        poisson_fft(int n_);
        ~poisson_fft();
        void init();
        void solve();
        void print();
        void init_mms();
        double l2_error_mms();
        void output(const char* filename,bool solution);
    private:
        double lambda(int j);
        fftw_plan plan_fwd;
        fftw_plan plan_bck;
};

#endif
