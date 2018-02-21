#include "sol_hammer_h.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Initializes the fourth-order Hammer--Hollingsworth solver, allocating
 * memory and setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
hammer_h::hammer_h(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k1b(new double[dof]),
    k2b(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
hammer_h::~hammer_h() {
    delete [] k2b;
    delete [] k1b;
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step with the fourth-order Hammer--Hollingsworth solver.
 * \param[in] dt the integration step. */
void hammer_h::step(double dt) {
    int iter=0;
    double delsq,d,*c;
    const double r36=sqrt(3)/6.;

    // Clear steps
    for(int i=0;i<dof;i++) k1[i]=k2[i]=0;

    do {

        // Check for too many iterations
        if(++iter>1000) {
            fputs("Too many iterations in IRK\n",stderr);
            exit(1);
        }

        // Perform update
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(0.25*k1[i]+(0.25-r36)*k2[i]);
        ff(t+dt*(0.5-r36),dq,k1b);
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*((0.25+r36)*k1[i]+0.25*k2[i]);
        ff(t+dt*(0.5+r36),dq,k2b);
        fcount+=2;

        // Find size of step from previous iteration
        delsq=0;
        for(int i=0;i<dof;i++) {
            d=k1[i]-k1b[i];delsq+=d*d;
            d=k2[i]-k2b[i];delsq+=d*d;
        }

	// Switch k1<->k1b and k2<->k2b array pointers. This will make k1 & k2
	// be used as the new values on the next iteration.
        c=k1b;k1b=k1;k1=c;
        c=k2b;k2b=k2;k2=c;
    } while(delsq>1e-25);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*0.5*(k1[i]+k2[i]);
    t+=dt;
}
