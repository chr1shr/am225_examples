#ifndef KDV_HH
#define KDV_HH

#include "sol_rk4.hh"

#include <fftw3.h>

/** This class has functions to specify the test Brusselator problem. */
class kdv : public rk4 {
    public:
        const int n;
        const int fftn;
        double* const der;
        double* const der3;
        kdv(int n_) : rk4(n_), n(n_), fftn(n/2+1),
            der(fftw_alloc_real(n)),
            der3(fftw_alloc_real(n)),
            c1(fftw_alloc_complex(fftn)),
            c2(fftw_alloc_complex(fftn)),
            plan1(fftw_plan_dft_r2c_1d(n,der,c1,FFTW_ESTIMATE)),
            plan2(fftw_plan_dft_c2r_1d(n,c1,der,FFTW_ESTIMATE)),
            plan3(fftw_plan_dft_c2r_1d(n,c2,der3,FFTW_ESTIMATE)) {}
        ~kdv() {
            fftw_destroy_plan(plan3);
            fftw_destroy_plan(plan2);
            fftw_destroy_plan(plan1);
            fftw_free(c2);
            fftw_free(c1);
            fftw_free(der3);
            fftw_free(der);
        }
        virtual void init();
        virtual void ff(double t_,double *in,double *out);
    private:
        fftw_complex* const c1;
        fftw_complex* const c2;
        fftw_plan plan1;
        fftw_plan plan2;
        fftw_plan plan3;
};

#endif
