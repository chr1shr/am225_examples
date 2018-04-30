#include "car_model.hh"

#include <cmath>

/** Evaluates the RHS for the discrete car model.
 * \param[in] t the dependent x variable in Hairer et al.'s notation.
 * \param[in] in the array containing y variable.
 * \param[in] out the function f(x,y). */
void car_model::ff(double t_,double *in,double *out) {
    for(int i=0;i<dof;i++)
        out[i]=1.-1./((i==dof-1?*in+L:in[i+1])-in[i]);
}

/** Sets up the initial conditions for the discrete car model. */
void car_model::init() {
    const double fac=L/dof,fac2=2*M_PI/dof;
    for(int i=0;i<dof;i++)
        q[i]=fac*i+14*exp(0.5*(-1-cos(fac2*i)));
}
