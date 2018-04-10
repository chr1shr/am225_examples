#include <cstdio>

#include "ho_transport.hh"

int main() {

    // Number of snapshots to output
    const int snaps=20;

    // Number of gridpoints
    const int m=256;

    // Integration timestep safety factor
    const double sf=0.2;

    // Type of integration. 0: Lax-Wendroff method, 1: Beam-Warming method
    const int type=3;

    // Create the transport simulation, and set up the initial condition
    ho_transport ho(m,1.0);
 //   ho.init_exp_sine();
    ho.init_step_function();

    // Integrate and save the solution snapshots to file
    ho.solve("ho3.out",snaps,1.0,sf,type);
}
