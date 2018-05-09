#include "sol_imp_onestep.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Initializes the one-step implicit method class, which incorporates several
 * different methods including the backward Euler method, and the implicit
 * midpoint method. The constructor allocates memory and sets constants.
 * \param[in] dof_ the number of degrees of freedom.
 * \param[in] alpha_ the fraction the f*/
imp_onestep::imp_onestep(int dof_,double alpha_) : sol_base(dof_), alpha(alpha_),
    dq(new double[dof]), k1(new double[dof]), k1b(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
imp_onestep::~imp_onestep() {
    delete [] k1b;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step with the one-step implicit method.
 * \param[in] dt the integration step.
 * \return True if the step was successful, false otherwise. */
bool imp_onestep::step(double dt) {
    int iter=0,done=0;
    double delsq,d,*c;

    // Clear steps
    for(int i=0;i<dof;i++) k1[i]=0;

    // Perform fixed point iteration. Loop until the termination criterion has
    // been reached several times.
    do {

        // Check for too many iterations
        if(++iter>1000) return false;

        // Perform update
        for(int i=0;i<dof;i++) dq[i]=q[i]+alpha*dt*k1[i];
        ff(t+alpha*dt,dq,k1b);
        fcount++;

        // Find size of step from previous iteration
        delsq=0;
        for(int i=0;i<dof;i++) {
            d=k1[i]-k1b[i];delsq+=d*d;
        }

        // Count if the termination criterion has been reached
        done=delsq<1e-25?done+1:0;

        // Switch k1<->k1b array pointers. This will make k1 be used as the new
        // values on the next iteration.
        c=k1b;k1b=k1;k1=c;
    } while(done<8);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*k1[i];
    t+=dt;
    return true;
}
