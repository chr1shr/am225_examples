#include "cubic_1d_fe.hh"

#include <cmath>

int main() {

#pragma omp parallel for schedule(dynamic) ordered
    for(int i=0;i<=30;i++) {
        int j=int(10*pow(100,i/30.));
        cubic_1d_fe cf(j);

        cf.init_mms();
        cf.solve();
#pragma omp ordered
        printf("%g %g\n",cf.h,cf.l2_norm_mms());
    }
}
