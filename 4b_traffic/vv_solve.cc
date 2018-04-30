#include "v_visco.hh"

#include <cmath>

int main() {

    // Case 1: high diffusion with epsilon=10^{-2}. Use a 1500 grid and simulate up to t=1.5,
    // outputting 15 snapshots
    {
        v_visco vv(1500,0.01);
        vv.initialize();
        vv.integrate(1.5,0.25,true,15);
        vv.output("vv_hi.out");
    }

    // Case 1: medium diffusion with epsilon=10^{-2.5}
    {
        v_visco vv(1500,sqrt(1e-5));
        vv.initialize();
        vv.integrate(1.5,0.25,true,15);
        vv.output("vv_md.out");
    }

    // Case 1: low diffusion with epsilon=10^{-3}
    {
        v_visco vv(1500,0.001);
        vv.initialize();
        vv.integrate(1.5,0.25,true,15);
        vv.output("vv_lo.out");
    }

}
