#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sol_rk4d.hh"
#include "osc.hh"

int main() {

    // Create the extrapolation integrator applied to the two-component
    // oscillator system. The argument sets the maximum j to be considered
    // in the T_{j,k} family of solutions.
    osc_extrap o(false,8);

    // Initialize the extrapolation integrator to use the harmonic sequence.
    // Other choices are the Romberg and Bulirsch sequences.
    o.init_harmonic();

    // Solve the system using a given T_{j,k} given by last two arguments.
    o.solve_fixed(8.,1000,true,1,1);
}
