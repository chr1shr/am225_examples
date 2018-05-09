#include "orbit.hh"

#include <cstdio>
#include <cmath>

/** Sets up the initial conditions for the ODE.
 * \param[in] q the array to write to. */
void orbit::orb_init(double *q) {
    const double r=a*(1+e);
    *q=0;
    q[1]=sqrt((1-e)/r);
    q[2]=r;
    q[3]=0;
    init_h=hamiltonian(q);
    //printf("# Initial H(p,q)=%.12g\n",init_h);
}

/** Prints the current state of the solution, plus the Hamiltonian.
 * \param[in] t_ the current simulation time.
 * \param[in] q the solution array. */
void orbit::orb_print(double t_,double *q) {
    printf("%g %g %g %g %g %.12g\n",t_,*q,q[1],q[2],q[3],hamiltonian(q)-init_h);
}

/** Computes the Hamiltonian.
 * \param[in] q the solution array.
 * \return The Hamiltonian. */
double orbit::hamiltonian(double *q) {
    return 0.5*(*q*(*q)+q[1]*q[1])-1/sqrt(q[2]*q[2]+q[3]*q[3]);
}
