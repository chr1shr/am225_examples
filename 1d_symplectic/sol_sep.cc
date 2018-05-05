#include "sol_sep.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Initializes the separable symplectic solver class. The constructor
 * allocates memory and sets constants.
 * \param[in] dof_ the number of degrees of freedom.
 * \param[in] mode the scheme to use. 0: first-order symplectic, 1: Ruth's
 *                 method, 2: fourth-order Ruth concatenation. */
sym_separable::sym_separable(int dof_,int mode) : sol_base(dof_),
    hdof(dof_>>1), steps(tab_steps[mode]), os(tab_offset[mode]),
    del(new double[hdof]) {}

/** The class destructor frees the dynamically allocated memory. */
sym_separable::~sym_separable() {
    delete [] del;
}

/** Performs an integration step.
 * \param[in] dt the integration step size.
 * \return True for successful completion. */
bool sym_separable::step(double dt) {
    double *cp=coeffs+os;fcount+=steps;

    for(int i=0;i<steps;i++) {

        // Update the momentum variables
        fp(q+hdof,del);
        for(int i=0;i<hdof;i++) q[i]+=*cp*dt*del[i];
        cp++;

        // Update the position variables
        fq(q,del);
        for(int i=0;i<hdof;i++) q[hdof+i]+=*cp*dt*del[i];
        cp++;
    }

    t+=dt;
    return true;
}

/** The number of steps for each different mode the class can run in. 0:
 * first-order symplectic, 1: Ruth's method, 2: fourth-order Ruth
 * concatenation. */
const int sym_separable::tab_steps[3]={1,3,6};

/** The offset in the coefficient table for the different modes. */
const int sym_separable::tab_offset[3]={0,2,8};

/** The coefficient table. */
double sym_separable::coeffs[20]={

    // First order
    1.,
    1.,

    // Third-order Ruth
    7./24,2./3,
    3./4,-2./3,
    -1./24,1,

    // Fourth-order Ruth concatenation
    7./48,1./3,
    3./8,-1./3,
    -1./48,1.,
    -1./48,-1./3,
    3./8,1./3,
    7./48,0.
};
