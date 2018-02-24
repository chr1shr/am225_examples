#include "schwarz.hh"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <limits>

/** Class for demonstrating the Schwarz algorithms on linear systems with
 * overlapping domains. It considers solving the Poisson equation on a domain
 * comprising of two overlapping square subdomains, each of which can be solved
 * quickly using the FFT.
 * \param[in] n_ the number of non-zero grid points in one direction on a
 *               square subdomain.
 * \param[in] (xo_,yo_) the offset (measure in grid points) of the second
 *                      square with respect to the first. */
schwarz::schwarz(int n_,int xo_,int yo_) : n(n_), xo(xo_), yo(yo_), ma(n+2+xo),
    na(n+2+yo), mna(ma*na), v(new double[mna]), f(new double[mna]),
    r(new double[mna]), grid1(n), grid2(n),
    tol(10*std::numeric_limits<double>::epsilon()/(grid1.h*grid1.h)),
    ntsq_a((2*n*n-(n-xo)*(n-yo))*tol*tol), ntsq_m(2*n*n*tol*tol) {
    double *p;

    // Set augmented solution and source terms grids to zero. This also helps
    // create clean output of the augmented grid, by setting unused parts to
    // zero.
    for(p=v;p<v+mna;p++) *p=0;
    for(p=f;p<f+mna;p++) *p=0;
}

/** The class destructor frees the dynamically allocated arrays. */
schwarz::~schwarz() {
    delete [] r;
    delete [] f;
    delete [] v;
}

/** Initializes the source term be a constant value of unity. */
void schwarz::init() {

    // Set up pointer to the first non-zero entry of the grid (skipping over
    // the boundary rows/columns)
    int i,j=0;
    double *fm=f+ma+1;

    // Initialize the rows of augmented grid that are only in the first
    // subdomain
    for(;j<yo;j++) for(i=0;i<n;i++) fm[j*ma+i]=1.;

    // Initialize the rows of augmented grid that ar common to both subdomains
    for(;j<n;j++) for(i=0;i<n+xo;i++) fm[j*ma+i]=1.;

    // Initialize the rows of augmented grid that are only in the first
    // subdomain
    for(;j<n+yo;j++) for(i=xo;i<n+xo;i++) fm[j*ma+i]=1.;
}

/** Solves the Poisson problem using the additive Schwarz method.
 * \param[in] verbose whether to print messages about convergence. */
void schwarz::solve_additive(bool verbose) {
    const double chr_factor=2;
    int i=0;

    // Pointers to square subdomains
    double *r1=r+ma+1,*r2=r+ma*(yo+1)+(xo+1),
           *v1=v+ma+1,*v2=v+ma*(yo+1)+(xo+1),

    // Do an initial computation of the residual
           rsq=compute_residual();
    if(verbose) printf("# Iter 0, residual %g\n",sqrt(rsq));

    while(rsq>ntsq_a) {

        // Check for too many iterations
        if(i++>1024) {
            fputs("Too many iterations\n",stderr);
            exit(1);
        }

        // Compute the solution updates on the two subdomains; this can be done
        // in parallel
#pragma omp parallel sections
        {
            {copy_to_square(r1,grid1.f);
            grid1.solve();}
#pragma omp section
            {copy_to_square(r2,grid2.f);
            grid2.solve();}
        }

        // Add the contributions
        for(int j=0;j<n-yo;j++) for(int k=0;k<n-xo;k++) v2[j*ma+k]*=chr_factor;
        add_from_square(grid1.v,v1);
        add_from_square(grid2.v,v2);
        for(int j=0;j<n-yo;j++) for(int k=0;k<n-xo;k++) v2[j*ma+k]*=1/chr_factor;

        rsq=compute_residual();
        if(verbose) printf("# Iter %d, residual %g\n",i,sqrt(rsq));
    }
}

/** Solves the Poisson problem using the multiplicative Schwarz method.
 * \param[in] verbose whether to print messages about convergence. */
void schwarz::solve_multiplicative(bool verbose) {

    // Pointers to square subdomains
    const int o1=ma+1,o2=ma*(yo+1)+(xo+1);
    int i=0;

    // Do an initial computation of the residual
    double rsq=residual_to_square(o1,grid1.f)
              +residual_to_square(o2,grid2.f);
    if(verbose) printf("# Iter 0, residual %g\n",sqrt(rsq));

    while(rsq>ntsq_m) {

        // Check for too many iterations
        if(i++>1024) {
            fputs("Too many iterations\n",stderr);
            exit(1);
        }

        // Compute the solution updates on the two subdomains. Unlike
        // the additive method, this cannot be done in parallel, since
        // the residual on the first subdomain is incorporated into the
        // second subdomain solve.
        grid1.solve();
        add_from_square(grid1.v,v+o1);
        rsq=residual_to_square(o2,grid2.f);
        grid2.solve();
        add_from_square(grid2.v,v+o2);

        rsq+=residual_to_square(o1,grid1.f);
        // Compute the rest of the residual, in preparation for the next
        // iteration
        if(verbose) printf("# Iter %d, residual %g\n",i,sqrt(rsq));
    }
}

/** Clears the solution array. */
void schwarz::clear_solution() {
    int i,j=0;
    double *vm=v+ma+1;
    for(;j<yo;j++) for(i=0;i<n;i++) vm[j*ma+i]=0;
    for(;j<n;j++) for(i=0;i<n+xo;i++) vm[j*ma+i]=0.;
    for(;j<n+yo;j++) for(i=xo;i<n+xo;i++) vm[j*ma+i]=0.;
}

/** Calculates the residual on a subdomain, and stores in the source term
 * of one of the subdomains.
 * \param[in] o the index in the augmented grid of the subdomain corner.
 * \param[in] g a pointer to the first entry in the subdomain grid.
 * \return The sum of square residuals. */
double schwarz::residual_to_square(int o,double *g) {
    double rsq=0;
    for(double *ge=g+n*n;g<ge;o+=ma-n) for(double *gr=g+n;g<gr;g++,o++) {
        *g=resid(o);
        rsq+=*g*(*g);
    }
    return rsq;
}

/** Computes the residual on the entire solution domain.
 * \return The sum of square residuals. */
double schwarz::compute_residual() {
    const double d=ma+1;

    // Set up pointer to the first non-zero entry of the grid (skipping over
    // the boundary rows/columns)
    int i,j=0;
    double rsq=0;

    // Loop over the grid points, compute the residual, and accumulate the L2 norm
    for(;j<yo;j++) for(i=0;i<n;i++) rsq+=store_resid(d+ma*j+i);
    for(;j<n;j++) for(i=0;i<n+xo;i++) rsq+=store_resid(d+ma*j+i);
    for(;j<n+yo;j++) for(i=xo;i<n+xo;i++) rsq+=store_resid(d+ma*j+i);
    return rsq;
}

/** Stores the residual at a grid point.
 * \param[in k the index of the grid point to consider.
 * \return The square of the residual. */
inline double schwarz::store_resid(int k) {
    r[k]=resid(k);
    return r[k]*r[k];
}

/** Calculates the residual at a gridpoint.
 * \return The residual. */
inline double schwarz::resid(int k) {
    return f[k]+(1/(grid1.h*grid1.h))*(v[-ma+k]+v[-1+k]-4*v[k]+v[1+k]+v[ma+k]);
}

/** Copies an n by n square of the residual vector in the augmented grid into
 * one of the subdomain grids.
 * \param[in] p a pointer to the first entry in the residual vector.
 * \param[in] g a pointer to the first entry in the subdomain grid. */
void schwarz::copy_to_square(double *p,double *g) {
    for(double *ge=g+n*n;g<ge;g+=n,p+=ma) memcpy(g,p,n*sizeof(double));
}

/** Add the solution from a subdomain grid onto the main solution vector in the
 * augmented grid.
 * \param[in] g a pointer to the first entry in the subdomain grid.
 * \param[in] p a pointer to the first entry in the residual vector. */
void schwarz::add_from_square(double *g,double *p) {
    for(double *ge=g+n*n;g<ge;p+=ma-n)
        for(double *gr=g+n;g<gr;) *(p++)+=*(g++);
}
