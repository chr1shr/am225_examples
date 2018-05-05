#include "orbit.hh"

#include <cstdio>
#include <cstring>
#include <cmath>

#include "omp.h"

int main() {

    // The semi-major axis of the orbit
    const double a=1.;

    // The eccentricity of the orbit
    const double e=1/3.;

    // The simulation duration, currently set to two complete orbits
    const double du=4*M_PI*a*sqrt(a);
    const bool exact_period=true;

    // An alternative simulation duration that is not an exact period
//    const double du=15.;
//    const bool exact_period=false;

    // Non-zero components of the initial condition
    double ref[4];

    // Compute reference solution
    if(exact_period) ref[2]=a*(1+e),ref[1]=sqrt((1-e)/ref[2]);
    else {
        orb_ruth4 sb(a,e);
        sb.solve_fixed(du,1000000);
        memcpy(ref,sb.q,4*sizeof(double));
    }

    // Allocate space for storing accuracy and function evaluations
    int fcount[707];
    double err[707],t0=omp_get_wtime();

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<=100;i++) {
        for(int j=0;j<7;j++) {

            // Compute the number of steps according to a power law. Minimum
            // steps of 100, and maximum steps of 500000.
            int ij=i+101*j,steps=int(100.*pow(5000.,0.01*i));

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
                double del,val=0.;
                for(int k=0;k<4;k++) del=sb->q[k]-ref[k],val+=del*del;
                err[ij]=sqrt(val);
                fcount[ij]=sb->fcount;
            } else fcount[ij]=-1;
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
