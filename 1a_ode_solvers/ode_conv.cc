#include "sol_euler.hh"
#include "sol_heun3.hh"
#include "sol_rk4.hh"
#include "sol_hammer_h.hh"
#include "brusselator.hh"

#include <cstdio>
#include <cmath>

#include "omp.h"

int main() {

    // Compute reference solution
    brus_rk4 br;
    br.solve_fixed(20.,200000);
    double ref0=br.q[0],ref1=br.q[1];

    // Allocate space for storing accuracy
    int fcount[505];
    double err[505],t0=omp_get_wtime();

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<=100;i++) {
        for(int j=0;j<5;j++) {

            // Compute the number of steps according to a power law. Minimum
            // steps of 100, and maximum steps of 100000.
            int steps=int(100.*pow(1000.,0.01*i));

            // Dynamically allocate a solver
            sol_base *sb;
            switch(j) {
                case 0: sb=new brus_euler();break;
                case 1: sb=new brus_ralston();break;
                case 2: sb=new brus_heun3();break;
                case 3: sb=new brus_rk4();break;
                default: sb=new brus_hammer_h();
            }

            // Perform integration and compute the difference to the reference
            // solution
            sb->solve_fixed(20.,steps);
            double dy0=ref0-sb->q[0],
                   dy1=ref1-sb->q[1];
            err[i+101*j]=sqrt(dy0*dy0+dy1*dy1);
            fcount[i+101*j]=sb->fcount;
            delete sb;
        }
    }

    // Output the results in a Gnuplot-readable format
    printf("# Time taken : %g ms\n",1e3*(omp_get_wtime()-t0));
    for(int j=0;j<5;j++) {
        for(int i=0;i<=100;i++) printf("%d %g\n",fcount[i+101*j],err[i+101*j]);
        puts("\n");
    }
}
