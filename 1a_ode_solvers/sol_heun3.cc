#include "sol_heun3.hh"

#include <cstdio>

/** Initializes the third-order Heun solver, allocating memory and setting
 * constants.
 * \param[in] dof_ the number of degrees of freedom. */
heun3::heun3(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
heun3::~heun3() {
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step with the third-order Heun solver.
 * \param[in] dt the integration step. */
void heun3::step(double dt) {

    // First RK step
    ff(t,q,k1);

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(1/3.)*dt*k1[i];
    ff(t+(1/3.)*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(2/3.)*dt*k2[i];
    ff(t+(2/3.)*dt,dq,k3);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*(0.25*k1[i]+0.75*k3[i]);
    t+=dt;fcount+=3;
}
