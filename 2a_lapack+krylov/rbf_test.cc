#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "rbf.hh"

int main() {

    rbf r(100,0);
    r.set_length_scale(0.5);

    r.init_random();
    r.solve_weights_lapack();
    r.output_points("rbf.pts");
    r.output_interpolant("rbf.fld",101,1.1);
}
