#include <cstdio>
#include <cmath>

#include <fftw3.h>
#include "omp.h"

const int n=24;
const double h=2*M_PI/n;
const double ninv=1./n;

int main() {
    int i;
    double re,im,x;

    // Allocate memory for FFTW input data, and convert to double precision
    double *fld=fftw_alloc_real(n),*der=fftw_alloc_real(n);
    for(i=0;i<n;i++) {
        x=i*h;
        fld[i]=exp(-cos(x));
    }

    // Allocate memory for complex FFTW output data
    int fftn=n/2+1;
    fftw_complex *c1=fftw_alloc_complex(fftn);

    // Make FFTW plans associated with performing this FFT
    fftw_plan plan1(fftw_plan_dft_r2c_1d(n,fld,c1,FFTW_ESTIMATE));
    fftw_plan plan2(fftw_plan_dft_c2r_1d(n,c1,der,FFTW_ESTIMATE));

    // Perform derivative
    fftw_execute(plan1);
    for(int i=0;i<fftn;i++) {
        re=c1[i][0];
        im=c1[i][1];
        c1[i][0]=-ninv*i*im;
        c1[i][1]=ninv*i*re;
    }
    fftw_execute(plan2);

    // Print solution and calculate L2 error
    double l2=0,eder,err;
    for(int i=0;i<n;i++) {
        x=i*h;

        // Exact derivative
        eder=sin(x)*exp(-cos(x));
        err=der[i]-eder;
        l2+=err*err;
        printf("%g %g %g %g %g\n",x,fld[i],der[i],eder,err);
    }
    printf("# L2 error: %g\n",sqrt(l2*ninv));

    // Remove memory that was dynamically allocated by FFTW
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(plan1);
    fftw_free(c1);
    fftw_free(der);
    fftw_free(fld);
}
