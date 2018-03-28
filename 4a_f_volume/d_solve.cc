#include <cstdio>

#include "diffuse.hh"

int main() {

    // Number of snapshots to output, and iterations between snapshot
    const int snaps=40,iters=20000;

    // Number of gridpoints
    const int m=512;

    // Integration timestep safety factor
    const double sf=0.2;

    // Type of integration. 0: finite-difference method, 1: finite-volume method
    const int type=0;

    // Create the diffusion simulation. Initialize the solution and the nu
    // table.
    diffuse d0(m);
    d0.init_step_function();
    d0.init_nu_array(0,sf);

    // Integrate and save the solution snapshots to file
    d0.solve("diff.out",snaps,iters,0);
}
