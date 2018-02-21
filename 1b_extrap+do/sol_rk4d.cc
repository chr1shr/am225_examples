#include "sol_rk4d.hh"

#include <cstdio>
#include <cstring>

/** Initializes the fourth-order Runge-Kutta solver, allocating memory and
 * setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
rk4d::rk4d(int dof_) : dof(dof_), fcount(0), t(0.), q(new double[dof]),
    dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k4(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
rk4d::~rk4d() {
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

/** Prints the current state of the solution. */
void rk4d::print(double t_,double *in) {
    printf("%g",t_);
    for(int i=0;i<dof;i++) printf(" %g",in[i]);
    puts("");
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void rk4d::solve_fixed(double duration,int steps,bool output,int d_steps) {

    // Set up initial condition and compute timestep. Use the t_den variable to
    // mark where the dense output has got to.
    init();
    double dt=duration/steps,t_den=0,dt_den;

    // Perform integration steps
    if(output) print(t,q);
    if(d_steps>0) {
        dt_den=duration/d_steps;
        print(t,q);
    }

    ff(t,q,k1);
    for(int i=0;i<steps;i++) {
        step(dt);

        // Do any dense output interpolation
        if(d_steps>0) {
            while(t_den+dt_den<t) {
                t_den+=dt_den;
                dense_output(1.+(t_den-t)/dt,dt);
                print(t_den,k3);
            }
        }

        // Move the required data for the next step into position
        memcpy(q,dq,dof*sizeof(double));
        memcpy(k1,k2,dof*sizeof(double));

        if(output) print(t,q);
    }
}

/** Computes a Hermite interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] theta the fraction of the timestep at which to evaluate the
 *                  interpolation.
 * \param[in] dt the length of the current timestep. */
void rk4d::dense_output(double theta,double dt) {
    double mth=1-theta;

    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k1, and the new derivative is in
    // k2
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+theta*dq[i]
             -theta*mth*((1-2*theta)*(dq[i]-q[i])+dt*(theta*k2[i]-mth*k1[i]));
}

/** Performs an integration step with the fourth-order Runge-Kutta solver.
 * \param[in] dt the integration step. */
void rk4d::step(double dt) {

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
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(1/6.)*(k1[i]+2*(k2[i]+k3[i])+k4[i]);

    // Reuse k2 to store the derivative at the new solution
    ff(t,dq,k2);
}
