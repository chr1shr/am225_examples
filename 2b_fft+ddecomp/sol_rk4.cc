#include "sol_rk4.hh"

#include <cstdio>

/** Initializes the fourth-order Runge-Kutta solver, allocating memory and
 * setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
rk4::rk4(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k4(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
rk4::~rk4() {
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step with the fourth-order Runge-Kutta solver.
 * \param[in] dt the integration step. */
void rk4::step(double dt) {

    // First RK step
    ff(t,q,k1);

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k1[i];
    ff(t+0.5*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+0.5*dt*k2[i];
    ff(t+0.5*dt,dq,k3);

    // Fourth RK step
    t+=dt;fcount+=4;
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k3[i];
    ff(t,dq,k4);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*(1/6.)*(k1[i]+2*(k2[i]+k3[i])+k4[i]);
}
