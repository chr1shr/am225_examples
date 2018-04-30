#include <cstdio>

#include "traffic.hh"

int main() {

    // Number of snapshots to output
    const int snaps=32;

    // Number of gridpoints
    const int m=256;

    // Integration timestep safety factor
    const double sf=0.1;

    // Create the traffic simulation, and set up the initial condition
    traffic tr(m,1.0);
 //   tr.init_exp_sine();
    tr.init_step_function();

    // Integrate and save the solution snapshots to file
    tr.solve("traf.out",snaps,2.0,sf);
}
