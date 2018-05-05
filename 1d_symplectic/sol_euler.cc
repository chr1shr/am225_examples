#include "sol_euler.hh"

#include <cstdio>

/** The class constructor when the number of degrees of freedom are known.
 * In this case the arrays for the
 * \param[in] dof_ the number of degrees of freedom. */
euler::euler(int dof_) : sol_base(dof_), k1(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
euler::~euler() {
    delete [] k1;
}

/** Performs an explicit Euler step.
 * \param[in] dt the integration step.
 * \return True for successful completion. */
bool euler::step(double dt) {
    ff(t,q,k1);
    for(int i=0;i<dof;i++) q[i]+=dt*k1[i];
    t+=dt;fcount++;
    return true;
}
