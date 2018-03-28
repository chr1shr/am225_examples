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
    const int type=1;

    // Create the diffusion simulation. Initialize the solution and the nu
    // table.
    diffuse ds(m);
    ds.init_step_function();
    ds.init_nu_array(type,sf);

    // Integrate and save the solution snapshots to file
    ds.solve("diff1.out",snaps,iters,type);
}
