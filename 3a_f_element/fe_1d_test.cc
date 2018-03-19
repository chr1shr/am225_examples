#include "cubic_1d_fe.hh"

int main() {

    cubic_1d_fe cf(32);

    cf.init_const();
    cf.g=0.;

    cf.solve();

    cf.print();
}
