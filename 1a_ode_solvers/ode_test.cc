#include "sol_euler.hh"
#include "sol_heun3.hh"
#include "sol_rk4.hh"
#include "sol_hammer_h.hh"
#include "brusselator.hh"

int main() {

    // Solve up to x=20 using the HH method with 500 steps
    brus_hammer_h bh;
    bh.solve_fixed(20.,500,true);
}
