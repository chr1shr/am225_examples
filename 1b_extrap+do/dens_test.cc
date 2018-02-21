#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sol_rk4d.hh"
#include "osc.hh"

/** The number of components in the ODE system. */
const int ns=2;

int main() {

    osc_rk4d o;
    o.solve_fixed(8.,80,false,750);
}
