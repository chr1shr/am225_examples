#include <cstdio>
#include <cmath>
#include <fftw3.h>
#include "omp.h"

// Maximum size of FFT to test
const int n_max=131072;

/** Measures the time taken for FFTW to perform a DFT on real data.
 * \param[in] (f,c) pointers to input and output arrays.
 * \param[in] n the number of points in the input array.
 * \param[in] flags the flags to use in the fftw_plan.
 * \return The wall clock time. */
double fftw_time(double *f,fftw_complex *c,int n,unsigned flags) {
    fftw_plan plan_dft(fftw_plan_dft_r2c_1d(n,f,c,flags));
    double t0=omp_get_wtime(),t1;int k=0;

    // Execute as many FFTs as will fit into one second
    do {
        fftw_execute(plan_dft);
        k++;t1=omp_get_wtime();
    } while(t1<t0+1);

    // Deallocate the plan and return the mean wall clock time
    fftw_destroy_plan(plan_dft);
    return (t1-t0)/k;
}

int main() {

    // Allocate some test data
    double *f=fftw_alloc_real(n_max);
    for(double *fp=f,x=0;fp<f+n_max;x+=2*M_PI/n_max)
        *(fp++)=sin(10.*pow(1.2,x))+cos(8*pow(0.8,x));

    // Allocate memory for complex FFTW output data
    fftw_complex *c=fftw_alloc_complex(n_max/2+1);

    // Time the FFTW execution using the three different planning routines 
    for(int n=16;n<=n_max;n<<=1)
        printf("%d %g %g %g %g\n",n,
               fftw_time(f,c,n,FFTW_ESTIMATE),
               fftw_time(f,c,n,FFTW_MEASURE),
               fftw_time(f,c,n,FFTW_PATIENT),
               fftw_time(f,c,n,FFTW_EXHAUSTIVE));

    // Free dynamically allocated memory
    fftw_free(c);
    fftw_free(f);
}
