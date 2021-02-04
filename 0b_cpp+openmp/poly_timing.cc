#include <cstdio>
#include <cmath>

#include "omp.h"

int main() {

    // The omp_get_wtime function returns the number of seconds of wall clock
    // time from an arbitrary datum
    double t0=omp_get_wtime(),t1,t2,t3,s=0,x;

    // Method 1 - direct multiplication
    for(x=0;x<1;x+=1e-8) s+=1+2*x+3*x*x+4*x*x*x+5*x*x*x*x;
    t1=omp_get_wtime();

    // Method 2 - power function in cmath library
    for(x=0;x<1;x+=1e-8) s+=1+2*x+3*pow(x,2)+4*pow(x,3)+5*pow(x,4);
    t2=omp_get_wtime();

    // Method 3 - Horner's method
    for(x=0;x<1;x+=1e-8) s+=1+x*(2+x*(3+x*(4+5*x)));

    // Print timing comparisons
    t3=omp_get_wtime();
    printf("Method 1 : %g s\n"
           "Method 2 : %g s\n"
           "Method 3 : %g s\n",t1-t0,t2-t1,t3-t2);
}
