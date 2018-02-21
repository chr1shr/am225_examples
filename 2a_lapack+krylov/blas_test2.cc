#include <cstdio>

// Include the BLAS functions
#include "blas.h"

#include "omp.h"

// Amount of time (in seconds) to run each test for
const double duration=2.;

int main() {

    // Maximum array size to try: set to 2^26, or 64 MB
    const int n=1<<26;

    // Create large arrays and fill with data
    double *a=new double[n],*b=new double[n],*c=new double[n];
    for(int i=0;i<n;i++) a[i]=b[i]=c[i]=1e-8*double(i);

    // Open output file
    FILE *f=fopen("btest2.out","w");
    if(f==NULL) {
        fputs("Can't open output file\n",stderr);
        return 1;
    }

    // Loop over a variety of different array sizes
    double alpha=0.12,beta=0.5,t0,t1,t2;
    char trans='n';
    int num_self,num_blas;
    for(int k=8;k<=(1<<13);k+=1+(k>>1)) {
        num_self=num_blas=0;

        // Test manual routine for doing b=b+alpha*a. First approach uses
        // a simple for-loop indexed on i.
        t0=omp_get_wtime();
        do {
            for(int i=0;i<k;i++) for(int j=0;j<k;j++)
               for(int e=0;e<k;e++) c[k*i+j]+=a[k*i+e]*b[k*e+j];
            t1=omp_get_wtime();num_self++;
        } while(t1<t0+duration);

        // Test corresponding BLAS routine. Due to C++/Fortran calling
        // conventions, all entries are passed as pointers, and there is an
        // underscore after the name.
        do {
            dgemm_(&trans,&trans,&k,&k,&k,&alpha,a,&k,b,&k,&beta,c,&k);
            t2=omp_get_wtime();num_blas++;
        } while(t2<t1+duration);

        // Print timing info, dividing test durations by the number of trials
        t0=(t1-t0)/num_self;
        t1=(t2-t1)/num_blas;
        printf("Array size %d: for(i) %g s, BLAS %g s\n",k,t0,t1);
        fprintf(f,"%d %g %g %d %d\n",k,t0,t1,num_self,num_blas);
    }

    // Free the dynamically allocated memory and close file
    fclose(f);
    delete [] a;
    delete [] b;
}
