#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "omp.h"

// Declare tolerance as a constant
const double tol=1e-14;

// Function to consider
double f(double x,double lam) {
    return lam-cos(x);
}

double ridders(double lam) {

    // Declare variables and their types explicitly
    double a=0,b=M_PI,fa=f(a,lam),fb=f(b,lam);
    double c,fc,d,fd,prev_d;

    int it=0;
    while(it<100) {

        // Midpoint evaluation
        c=0.5*(b+a);
        fc=f(c,lam);

        // Ridders calculation
        d=c-(c-a)*(fa>fb?-fc:fc)/sqrt(fc*fc-fa*fb);
        fd=f(d,lam);

        // Check for acceptable solution
        if(it>0&&fabs(prev_d-d)<tol) return d;

        // Update counters, using C++ 'post-increment' operator on the it
        // variable
        it++;
        prev_d=d;

        // Reduce interval size
        if(fd<0) {
            a=d;fa=fd;
            if(fc>=0) {b=c;fb=fc;}
        } else {
            b=d;fb=fd;
            if(fc<0) {a=c;fa=fc;}
        }
    }
    fputs("# Too many iterations\n",stderr);
    exit(1);
}

int main() {

    // Total table size, and lambda increment
    const int ts=1000000;
    const double dlam=1.98/(ts-1);

    // Initialize table
    double xv[ts],t0=omp_get_wtime(),dt;

    // Time the table construction
    for(int l=0;l<ts;l++) xv[l]=ridders(-0.99+dlam*l);

    // Print timing results
    dt=omp_get_wtime()-t0;
    printf("Time: %.4g s (total)\nTime: %g microseconds (per value)\n",dt,1e6*dt/ts);

    // Print sum of terms. This just removes the compiler warning about the xv
    // values being unused.
    double s=0;
    for(int l=0;l<ts;l++) s+=xv[l];
    printf("Sum: %g\n",s);
}
