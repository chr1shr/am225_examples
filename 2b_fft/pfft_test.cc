#include "poisson_fft.hh"

int main() {
    poisson_fft pf(32);

    pf.init_mms();
    pf.solve();

    pf.output("pfft.sol",true);
    pf.output("pfft.inp",false);
}
