#include <cstdio>

#include "poisson_fft.hh"

int main() {

    for(int n=10;n<=1024;n+=n/3) {
        poisson_fft pf(n);
        pf.init_mms();
        pf.solve();
        printf("%d %g\n",n,pf.l2_error_mms());
    }
}
