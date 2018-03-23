#include "conj_grad.hh"

#include <cstdio>
#include <cmath>
#include <limits>

/** Reasonable guess for the tolerance criterion in the conjugate gradient
 * solver, assuming solution vector is order unity. */
const double conj_grad_tol=std::numeric_limits<double>::epsilon()*1e3;
const double conj_grad_tol_sq=conj_grad_tol*conj_grad_tol;

/** Initializes the conjugate gradient solver class, setting constants and
 * allocating memory.
 * \param[in] dof_ the degrees of freedom. */
conj_grad::conj_grad(int dof_) : dof(dof_), one(1), x(new double[dof]),
    b(new double[dof]), rk(new double[dof]), pk(new double[dof]),
    zk(new double[dof]), yk(NULL) {}

/** The class destructor frees the dynamically allocated memory. */
conj_grad::~conj_grad() {
    if(yk!=NULL) delete [] yk;
    delete [] zk;
    delete [] pk;
    delete [] rk;
    delete [] b;
    delete [] x;
}

/** Solves the linear system using the standard conjugate gradient method.
 * \param[in] verbose whether to output messages about residual. */
void conj_grad::solve(bool verbose) {
    int i,k=0;

    // Setup for conjugate gradient algorithm
    for(i=0;i<dof;i++) x[i]=0.;
    copy(b,rk);
    copy(b,pk);
    double rsq=dot(rk,rk),nrsq,nu;

    // Conjugate gradient algorithm, using BLAS commands for high efficiency
    while(rsq>conj_grad_tol_sq*dof) {
        if(verbose) printf("# Iter %d, residual %g\n",k++,sqrt(rsq));
        mul_A(pk,zk);
        nu=rsq/dot(pk,zk);
        daxpy(nu,pk,x);
        daxpy(-nu,zk,rk);
        dscal((nrsq=dot(rk,rk))/rsq,pk);
        daxpy(1.,rk,pk);
        rsq=nrsq;
    }
    if(verbose) printf("# Iter %d, residual %g\n",k,sqrt(rsq));
}

/** Solves the linear system using the standard conjugate gradient method. */
void conj_grad::solve_pre(bool verbose) {
    int i,k=0;

    // Allocate extra array of not already available
    if(yk==NULL) yk=new double[dof];

    // Setup for conjugate gradient algorithm
    for(i=0;i<dof;i++) x[i]=0.;
    copy(b,rk);
    M_inv(b,pk);
    copy(pk,yk);
    double rsq=dot(yk,rk),nrsq,nu;

    // Conjugate gradient algorithm, using BLAS commands for high efficiency
    while(dot(rk,rk)>conj_grad_tol_sq*dof) {
        if(verbose) printf("# Iter %d, residual %g\n",k++,sqrt(rsq));
        mul_A(pk,zk);
        nu=rsq/dot(pk,zk);
        daxpy(nu,pk,x);
        daxpy(-nu,zk,rk);
        M_inv(rk,yk);
        dscal((nrsq=dot(yk,rk))/rsq,pk);
        daxpy(1.,yk,pk);
        rsq=nrsq;
    }
    if(verbose) printf("# Iter %d, residual %g\n",k,sqrt(rsq));
}

/** This default preconditioning routine just applies the identity matrix.
 * \param[in] in the input vector.
 * \param[in] out the output vector after applying the preconditioning M^{-1},
 *                which in this case is the identity. */
void conj_grad::M_inv(double *in,double *out) {
    copy(in,out);
}
