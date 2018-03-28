#include <cstdio>

#include "diffuse.hh"

int main() {

    // Number of snapshots to output, and iterations between snapshot
    const int snaps=40,iters=20000;

    // Number of gridpoints
    const int m=512;

    // Integration timestep safety factor
    const double sf=0.2;

    // Create two diffusion simulations, to run with the two update methods. Initialize
    // both of them;
    diffuse d0(m),d1(m);
    d0.init_step_function();
    d1.init_step_function();
    d0.init_nu_array(0,sf);
    d1.init_nu_array(1,sf);

    // Integrate both classes
    double tsnap=d0.dt*iters;
    printf("0 %g %g\n",d0.integral(),d1.integral());
    for(int i=1;i<=snaps;i++) {

        // Step both simulations forward
#pragma omp parallel sections
        {
#pragma omp section
            for(int k=0;k<iters;k++) d0.step_forward();
#pragma omp section
            for(int k=0;k<iters;k++) d1.step_forward();
        }

        // Compute the integral of both solutions
        printf("%g %g %g\n",iters*d0.dt,d0.integral(),d1.integral());
    }
}
