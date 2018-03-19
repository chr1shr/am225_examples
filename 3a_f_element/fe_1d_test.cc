#include <cmath>

#include "cubic_1d_fe.hh"

int main() {

    cubic_1d_fe cf(90);
    cf.g=0.5;//exp(-1)*5*M_PI;

    cf.init_const();

//    cf.print_matrix();
    cf.solve();

    cf.print();
}
