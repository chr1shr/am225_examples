#include <cstdio>
#include <cmath>

#include <fftw3.h>
#include "omp.h"

const int n=8192;

int main() {
    int i;
    float e[n];
    double re,im;

    // Read in the binary data in single precision
    FILE *fp=fopen("bob.raw","rb");
    fread(e,sizeof(float),n,fp);
    fclose(fp);

    // Allocate memory for FFTW input data, and convert to double precision
    double *fld=(double*) fftw_malloc(sizeof(double)*n);
    for(i=0;i<n;i++) {
        fld[i]=e[i];
    }

    // Allocate memory for complex FFTW output data
    int fftn=n/2+1;
    fftw_complex *c1=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fftn);

    // Make FFTW plan associated with performing this FFT
    fftw_plan plan_dft(fftw_plan_dft_r2c_1d(n,fld,c1,FFTW_ESTIMATE));
    fftw_execute(plan_dft);

    // Optional timing routine
    /*double t0=omp_get_wtime();
    for(int i=0;i<50000;i++) fftw_execute(plan_dft);
    printf("Time %g ms\n",1e3/50000*(omp_get_wtime()-t0));
    return 1;*/

    // Output magnitudes of each term
    for(i=0;i<fftn;i++) {
        re=c1[i][0];
        im=c1[i][1];
        printf("%g %g\n",1/(4096/44000.)*i,sqrt(re*re+im*im));
    }
}    
