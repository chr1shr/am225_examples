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
    d0.init_exp_sine();
    d1.init_exp_sine();
    d0.init_nu_array(0,sf);
    d1.init_nu_array(1,sf);

    // Integrate both classes
    double tsnap=d0.dt*iters,id0=d0.integral(),id1=d1.integral();
    printf("0 %g %g 0 0\n",id0,id1);
    for(int i=1;i<=snaps;i++) {

        // Step both simulations forward
#pragma omp parallel sections
        {
#pragma omp section
            for(int k=0;k<iters;k++) d0.step_forward<0>();
#pragma omp section
            for(int k=0;k<iters;k++) d1.step_forward<1>();
        }

        // Compute the integral of both solutions
        double cd0=d0.integral(),cd1=d1.integral();
        printf("%g %g %g %g %g\n",i*tsnap,cd0,cd1,cd0-id0,cd1-id1);
    }
}
