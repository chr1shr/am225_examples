#include <cstdio>

#ifdef _OPENMP
#include "omp.h"
#endif

void dummy_computation(int *p,int l) {
    printf("Computing on address %p\n",(void*) p);
    for(int k=0;k<1000000000;k++) {
        *p<<=1;*p^=k;*p+=l;
    }
}

int main() {
    int t=omp_get_max_threads(),
        *s=new int[t];

    for(int k=0;k<t;k++) s[k]=k;
    double t0=omp_get_wtime(),t1;

#pragma omp parallel for
    for(int k=0;k<t;k++) dummy_computation(s+k,k);
    t1=omp_get_wtime();

#pragma omp parallel
    {
        int q[1];
        dummy_computation(q,omp_get_thread_num());
    }

    printf("# Time for method 1 : %g s\n"
           "# Time for method 2 : %g s\n",t1-t0,omp_get_wtime()-t0);
}
