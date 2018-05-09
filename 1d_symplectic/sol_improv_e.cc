#include "sol_improv_e.hh"

#include <cstdio>

/** Initializes the second-order improved Euler solver, allocating memory and
 * setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
improv_e::improv_e(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
improv_e::~improv_e() {
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step of the second-order Ralston method.
 * \param[in] dt the integration step.
 * \return True for successful completion. */
bool improv_e::step(double dt) {

    // First RK step
    ff(t,q,k1);

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k1[i];
    ff(t+dt,dq,k2);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=0.5*dt*(k1[i]+k2[i]);
    t+=dt;fcount+=2;

    return true;
}
