#include "sol_ralston.hh"

#include <cstdio>

/** Initializes the second-order Ralston solver, allocating memory and setting
 * constants.
 * \param[in] dof_ the number of degrees of freedom. */
ralston::ralston(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
ralston::~ralston() {
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step of the second-order Ralston method.
 * \param[in] dt the integration step. */
void ralston::step(double dt) {

    // First RK step
    ff(t,q,k1);

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(2/3.)*dt*k1[i];
    ff(t+(2/3.)*dt,dq,k2);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*(0.25*k1[i]+0.75*k2[i]);
    t+=dt;fcount+=2;
}
