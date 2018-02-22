#include "sol_base.hh"

#include <cstdio>

/** The class constructor initializes the class constants and dynamically
 * allocates memory for the solution.
 * \param[in] dof_ the number of degrees of freedom. */
sol_base::sol_base(int dof_) : dof(dof_), fcount(0), t(0.), q(new double[dof]) {}

/** The class constructor frees the dynamically allocated memory. */
sol_base::~sol_base() {
    delete [] q;
}

/** Prints the current state of the solution. */
void sol_base::print() {
    printf("%g",t);
    for(int i=0;i<dof;i++) printf(" %g",q[i]);
    puts("");
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void sol_base::solve_fixed(double duration,int steps,bool output) {

    // Set up initial condition and compute timestep
    init();
    double dt=duration/steps;

    // Perform integration steps
    if(output) print();
    for(int i=0;i<steps;i++) {
        step(dt);
        if(output) print();
    }
}
