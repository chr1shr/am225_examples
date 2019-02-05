#include <cstdio>

#ifdef _OPENMP
#include "omp.h"
#endif

/** Performs a complicated calculation on a section of memory containing four
 * integers.
 * \param[in] a pointer to the memory. */
void dummy_computation(int *p) {
    int l,t=omp_get_thread_num();
    printf("Thread %d computing on address %p\n",t,(void*) p);
    for(int k=0;k<10000000;k++) {
        l=*p&3;p[l]^=k+t;p[3-l]<<=1;
    }
}

int main() {

    // Create a small array and fill it with initial data
    int n=omp_get_max_threads(),*s=new int[4*n];
    for(int k=0;k<4*n;k++) s[k]=k;

    // Time how long it takes for multiple threads to perform the calculation
    // on nearby memory, potentially within a single cache line
    double t0=omp_get_wtime(),t1;
    puts("Method A - operating on nearby memory:");
#pragma omp parallel for
    for(int k=0;k<n;k++) dummy_computation(s+4*k);
    t1=omp_get_wtime();
    delete [] s;

    // Time how long it takes for multiple thread to perform the calculation if
    // they allocate their own memory, which is guaranteed to be in separate
    // locations
    puts("\nMethod B - operating on separate memory:");
#pragma omp parallel
    {
        int *q=new int[4];
        dummy_computation(q);
        delete [] q;
    }

    // Report the times for the two computations
    printf("\nTime for method A : %g s\n"
           "Time for method B : %g s\n",t1-t0,omp_get_wtime()-t1);
}
