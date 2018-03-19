#include <cmath>

#include "cubic_1d_fe.hh"

int main() {

    cubic_1d_fe cf(32);
    cf.g=-0.05;

    cf.init_slope();

//    cf.print_matrix();
    cf.solve();

    cf.print();
}
