#include "poisson_fft.hh"

int main() {
    poisson_fft pf(32);

    pf.init();
    pf.solve();

    pf.output_solution("pfft.sol");
    pf.output_source("pfft.src");
}
