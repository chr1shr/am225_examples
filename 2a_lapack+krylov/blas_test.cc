#include <cstdio>

// Include the BLAS functions
#include "blas.h"

#include "omp.h"

// Amount of time (in seconds) to run each test for
const double duration=4.;

int main() {

    // Maximum array size to try: set to 2^26, or 64 MB
    const int n=1<<26;

    // Create large arrays and fill with data
    double *a=new double[n],*b=new double[n];
    for(int i=0;i<n;i++) a[i]=b[i]=double(i);

    // Open output file
    FILE *f=fopen("btest.out","w");
    if(f==NULL) {
        fputs("Can't open output file\n",stderr);
        return 1;
    }

    // Loop over a variety of different array sizes
    double alpha=0.12,t0,t1,t2,t3;
    int l=1,num_for1,num_for2,num_blas;
    for(int k=4096;k<=n;k<<=1) {
        num_for1=num_for2=num_blas=0;

        // Test manual routine for doing b=b+alpha*a. First approach uses
        // a simple for-loop indexed on i.
        t0=omp_get_wtime();
        do {
            for(int i=0;i<k;i++) b[i]+=alpha*a[i];
            t1=omp_get_wtime();num_for1++;
        } while(t1<t0+duration);

        // Second approach uses pointers
        do {
            for(double *ap=a,*bp=b;ap<a+k;ap++,bp++) *bp+=*ap*alpha;
            t2=omp_get_wtime();num_for2++;
        } while(t2<t1+duration);

        // Test corresponding BLAS routine. Due to C++/Fortran calling
        // conventions, all entries are passed as pointers, and there is an
        // underscore after the name.
        do {
            daxpy_(&k,&alpha,a,&l,b,&l);
            t3=omp_get_wtime();num_blas++;
        } while(t3<t2+duration);

        // Print timing info, dividing test durations by the number of trials
        t0=(t1-t0)/num_for1;
        t1=(t2-t1)/num_for2;
        t2=(t3-t2)/num_blas;
        printf("Array size %d: for(i) %g s, for(ptr) %g s, BLAS %g s\n",k,t0,t1,t2);
        fprintf(f,"%d %g %g %g %d %d %d\n",k,t0,t1,t2,num_for1,num_for2,num_blas);
    }

    // Free the dynamically allocated memory and close file
    fclose(f);
    delete [] a;
    delete [] b;
}
