#include <cstdio>

#include "poisson_fft.hh"

int main() {

    // Loop over a range of grid sizes. Note that while FFTW achieves better
    // performance when the grid size is a power of two, it will achieve
    // n*log(n) scaling for all n.
    for(int n=10;n<=2048;n+=n>>1) {
        poisson_fft pf(n);

        // Initialize the source term for the manufactured solution
        pf.init_mms();

        // Solve the system and compute the L2 error to the analytical
        // manufactured solution
        pf.solve();
        printf("%d %g\n",n,pf.l2_error_mms());
    }
}
