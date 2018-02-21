#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sol_extrap.hh"
#include "osc.hh"

/** The final integration time. */
const double tend=5;

int main() {

    // Loop over a variety of T_{j,k} possibilities
    for(int k=1;k<=7;k++) {
        for(int j=k;j<=7;j++) {

            // Allocate a solver and set up the Bulirsch sequence
            osc_extrap o(true,10);
            o.init_bulirsch();

            // Perform integration and compute the difference to the reference
            // solution
            o.solve_fixed(tend,100,false,j,k);

            // Calculate convergence to exact solution
            double dy0=o.sol0(tend)-o.q[0],
                   dy1=o.sol1(tend)-o.q[1];
            printf("%d %d %d %g\n",j,k,o.fcount,sqrt(dy0*dy0+dy1*dy1));
        }
        puts("\n");
    }
}
