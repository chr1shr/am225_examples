#include "poisson_fft.hh"

int main() {
    poisson_fft pf(64);

    pf.init();
    pf.solve();

    pf.output_solution("pfft.sol");
    pf.output_source("pfft.src");
}
