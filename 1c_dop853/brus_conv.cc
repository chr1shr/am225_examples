#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "dop853.hh"
#include "functions.hh"

/** The number of components in the ODE system. */
const int ns=2;

/** The initial integration time. */
const double tstart=0;

/** The final integration time. */
const double tend=20;

int main() {

    // Construct the test class, set the level of accuracy, and carry out
    // the reference ODE solution
    brusselator o;
    o.atoli=1e-15;
    o.rtoli=1e-15;
    o.solve(tstart,tend,0,NULL);
    double ref0=o.w[0],ref1=o.w[1];

    // Initialize arrays for stor convergence info
    int fcount[101];
    double err[101];

#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<=100;i++) {

        // Dynamically allocate a solver
        brusselator br;
        br.atoli=1e-2*pow(0.1,0.1*i);
        br.rtoli=1e-2*pow(0.1,0.1*i);

        // Perform integration and compute the difference to the reference
        // solution
        br.solve(0,20,0,NULL);
        double dy0=ref0-br.w[0],
               dy1=ref1-br.w[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=br.nff;
    }

    // Output the results in a Gnuplot-readable format
    for(int i=0;i<=100;i++) printf("%d %g\n",fcount[i],err[i]);
}
