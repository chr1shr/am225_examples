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
    hdof(dof_>>1), steps(tab_steps[mode]), cstart(coeffs+tab_offset[mode]),
    del(new double[hdof]) {}

/** The class destructor frees the dynamically allocated memory. */
sym_separable::~sym_separable() {
    delete [] del;
}

/** Performs an integration step.
 * \param[in] dt the integration step size.
 * \return True for successful completion. */
bool sym_separable::step(double dt) {

    for(double *cp=cstart,*ce=cp+2*steps;cp<ce;) {

        // Update the momentum variables. For the concatenated method, assume
        // fp does not need to be called, apart from on the very first step.
        if(steps<5||cp>cstart||fcount==0) fp(q+hdof,del);
        for(int i=0;i<hdof;i++) q[i]+=*cp*dt*del[i];
        cp++;

        // Update the position variables
        fq(q,del);
        for(int i=0;i<hdof;i++) q[hdof+i]+=*cp*dt*del[i];
        cp++;
    }

    // Deal separately with the final step of the concatenated method
    if(steps==5) {
        if(fcount==0) fcount++;
        fp(q+hdof,del);
        for(int i=0;i<hdof;i++) q[i]+=7./48*dt*del[i];
    }

    t+=dt;fcount+=2*steps;
    return true;
}

/** The number of steps for each different mode the class can run in. 0:
 * first-order symplectic, 1: Ruth's method, 2: fourth-order Ruth
 * concatenation. */
const int sym_separable::tab_steps[3]={1,3,5};

/** The offset in the coefficient table for the different modes. */
const int sym_separable::tab_offset[3]={0,2,8};

/** The coefficient table. */
double sym_separable::coeffs[18]={

    // First order
    1.,
    1.,

    // Third-order Ruth
    7./24,2./3,
    3./4,-2./3,
    -1./24,1,

    // Fourth-order Ruth concatenation. (The sixth step, which only has one
    // non-zero term, is dealt with separately.)
    7./48,1./3,
    3./8,-1./3,
    -1./48,1.,
    -1./48,-1./3,
    3./8,1./3,
};
