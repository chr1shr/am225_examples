#include "poisson.hh"

#include <cstdio>

/** The class constructor sets required constants.
 * \param[in] dof the number of grid points. */
poisson::poisson(int dof_) : conj_grad(dof_), dx(1./(dof-1)),
    xsp2(1/(dx*dx)) {}

/** Initializes the source term for the matrix problem. */
void poisson::init() {
    *b=b[dof-1]=0;
    for(int i=1;i<dof-1;i++) b[i]=1.;
}

/** Prints the solution vector in text format. */
void poisson::print_solution() {
    for(int i=0;i<dof;i++) printf("%g %g\n",i*dx,x[i]);
}

/** Performs the matrix-vector operation required for the Krylov solver.
 * \param[in] in the input vector.
 * \param[in] out the output vector, after multiplying by the matrix. */
void poisson::mul_A(double *in,double *out) {

    // Dirichlet condition on left edge
    *out=-2*xsp2*(*in);

    // Second order centered finite-difference
    for(int i=1;i<dof-1;i++)
        out[i]=xsp2*(in[i-1]-2*in[i]+in[i+1]);

    // Dirichlet condition on right edge
    out[dof-1]=-2*xsp2*in[dof-1];
}
