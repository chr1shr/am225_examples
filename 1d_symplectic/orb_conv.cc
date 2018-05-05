#include "orbit.hh"

#include <cstdio>
#include <cmath>

#include "omp.h"

int main() {

    // The semi-major axis of the orbit
    const double a=1.;

    // The eccentricity of the orbit
    const double e=1/3.;

    // The simulation duration, currently set to one complete orbit
    const double du=6*M_PI*a*sqrt(a);

    // Non-zero components of the initial condition
    const double ref2=a*(1+e),ref1=sqrt((1-e)/ref2);

    // Allocate space for storing accuracy
    int fcount[707];
    double err[707],t0=omp_get_wtime();

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<=100;i++) {
        for(int j=0;j<7;j++) {

            // Compute the number of steps according to a power law. Minimum
            // steps of 100, and maximum steps of 500000.
            int steps=int(100.*pow(5000.,0.01*i));

            // Dynamically allocate a solver
            sol_base *sb;
            switch(j) {
                case 0: sb=new orb_euler(a,e);break;
                case 1: sb=new orb_back(a,e);break;
                case 2: sb=new orb_improv_e(a,e);break;
                case 3: sb=new orb_imp_mid(a,e);break;
                case 4: sb=new orb_fo_symplectic(a,e);break;
                case 5: sb=new orb_ruth3(a,e);break;
                default: sb=new orb_ruth4(a,e);
            }

            // Perform integration and compute the difference to the reference
            // solution, as well as the integration steps taken. In the case
            // when the integration terminated early, mark it as counting
            // negative number of steps to indicate an invalid result.
            if(sb->solve_fixed(du,steps)) {
                double d1=ref1-sb->q[1],d2=ref2-sb->q[2];
                err[i+101*j]=sqrt(d1*d1+d2*d2+sb->q[0]*sb->q[0]+sb->q[3]*sb->q[3]);
                fcount[i+101*j]=sb->fcount;
            } else fcount[i+101*j]=-1;
            delete sb;
        }
    }

    // Output the results in a Gnuplot-readable format
    printf("# Time taken : %g ms\n",1e3*(omp_get_wtime()-t0));
    for(int j=0;j<7;j++) {
        for(int i=0;i<=100;i++) {
            int ij=i+101*j;
            if(fcount[ij]>0) printf("%d %g\n",fcount[ij],err[ij]);
            else puts("NaN NaN");
        }
        puts("\n");
    }
}
