#include "cubic_1d_fe.hh"

int main() {

#pragma omp parallel for schedule(dynamic)
    for(int i=10;i<1000;i+=i>>1) {
        cubic_1d_fe cf(i);

        cf.init_mms();
        cf.solve();
#pragma omp ordered
        printf("%g %g\n",cf.h,cf.l2_norm_mms());
    }
}
