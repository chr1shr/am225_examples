#include "orbit.hh"

#include <cstdio>
#include <cmath>

/** Evaluates the function f(x,y) on the RHS of the ODE.
 * \param[in] t the dependent x variable in Hairer et al.'s notation.
 * \param[in] in the array containing y variable.
 * \param[in] out the function f(x,y). */
void orbit::orb_ff(double t_,double *in,double *out) {
    double fac=1/(in[2]*in[2]+in[3]*in[3]);
    fac*=sqrt(fac);
    *out=-in[2]*fac;
    out[1]=-in[3]*fac;
    out[2]=*in;
    out[3]=in[1];
}

/** Sets up the initial conditions for the ODE.
 * \param[in] q the array to write to. */
void orbit::orb_init(double *q) {
    *q=0.;
    q[1]=0.7;
    q[2]=1.;
    q[3]=0.;
    init_h=hamiltonian(q);
    printf("# Initial H(p,q)=%.12g\n",init_h);
}

/** Prints the current state of the solution, plus the Hamiltonian.
 * \param[in] t_ the current simulation time.
 * \param[in] q the solution array. */
void orbit::orb_print(double t_,double *q) {
    printf("%g %g %g %g %g %.10g\n",t_,*q,q[1],q[2],q[3],hamiltonian(q)-init_h);
}

/** Computes the Hamiltonian.
 * \param[in] q the solution array.
 * \return The Hamiltonian. */
double orbit::hamiltonian(double *q) {
    return 0.5*(*q*(*q)+q[1]*q[1])-1/sqrt(q[2]*q[2]+q[3]*q[3]);
}
