#include "sol_euler.hh"
#include "sol_improv_e.hh"
#include "sol_imp_onestep.hh"
#include "orbit.hh"

int main() {

    // Solve up to x=20
    orb_improv_e bh;
    bh.solve_fixed(10.,800,true);
}
