#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sol_extrap.hh"
#include "osc.hh"

/** The number of components in the ODE system. */
const int ns=2;

/** The initial integration time. */
const double tstart=0;

/** The final integration time. */
const double tend=5;

int main() {

    // Initialize arrays for storing convergence info
    int fcount[101];
    double err[101];

#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<=100;i++) {

        // Allocate a solver
        osc_extrap o(false,10);
        o.init_bulirsch();

        // Perform integration and compute the difference to the reference
        // solution
        int steps=int(100.*pow(1000.,i*0.01));
        o.solve_fixed(tend,steps,false,6,6);

        double dy0=o.sol0(tend)-o.q[0],
               dy1=o.sol1(tend)-o.q[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=o.fcount;
    }

    // Output the results in a Gnuplot-readable format
    for(int i=0;i<=100;i++) printf("%d %g\n",fcount[i],err[i]);
}
