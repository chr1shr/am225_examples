#include <cstdio>
#include <cmath>
#include <fftw3.h>
#include "omp.h"

const int n=8192;

int main() {

    // Read in the binary data in single precision
    float e[n];
    FILE *fp=fopen("bob.raw","rb");
    fread(e,sizeof(float),n,fp);
    fclose(fp);

    // Allocate memory for FFTW input data, and convert sound
    // sample to double precision
    double *f=fftw_alloc_real(n),re,im;
    for(int i=0;i<n;i++) f[i]=e[i];

    // Allocate memory for complex FFTW output data
    int fftn=n/2+1;
    fftw_complex *c=fftw_alloc_complex(n);

    // Make FFTW plan, and perform the transform
    fftw_plan plan_dft(fftw_plan_dft_r2c_1d(n,f,c,FFTW_ESTIMATE));
    fftw_execute(plan_dft);

    // Output magnitudes of each term
    for(int i=0;i<fftn;i++) {
        re=c[i][0];
        im=c[i][1];
        printf("%g %g\n",44000./(n/2)*i,sqrt(re*re+im*im));
    }

    // Free dynamically allocated memory
    fftw_destroy_plan(plan_dft);
    fftw_free(c);
    fftw_free(f);
}
