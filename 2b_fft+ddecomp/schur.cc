#include <cstdio>
#include "schur.hh"

/** Construct given the width of each subdomain.
 * \param[in] n_ : The width of each subdomain in units. */
schur::schur(int n_) : conj_grad(n_-1), n(n_-1), nn(n*n), ng(n), h(1./n_), ih2(n_*n_),
    grid1(n), grid2(n), f1(new double[nn]), f2(new double[nn]), fg(new double[ng])
{}

schur::~schur()
{
    delete [] fg;
    delete [] f2;
    delete [] f1;
}

/** Solve Poisson using the Schur complement method.
 * \param[in] f : A function handle to the source term. */
void schur::solve(const std::function<double(double,double)>& f)
{
    // Evaluate the source term at the gridpoints
    for (int i=0; i<n; ++i) {
        int ii = n-i;
        fg[i] = f(0, ii*h);
        for (int j=0; j<n; ++j) {
            f1[n*i+j] = grid1.f[n*i+j] = f(-1+(j+1)*h, ii*h);
            f2[n*i+j] = grid2.f[n*i+j] = f(   (j+1)*h, ii*h);
        }
    }

    // First, solve the Schur complement system for the unknowns on the glue
    //
    //    S*ug = fg - Ag1*A11^-1*f1 - Ag2*A22^-1*f2
    //
    // where S = Agg - Ag1*A11^-1*A1g - Ag2*A22^-1*A2g.
    // Here Aii^-1 means "call the Poisson FFT solver on subdomain i."

    grid1.solve(); // Compute A11^-1*f1
    grid2.solve(); // Compute A22^-1*f2

    // Compute b = fg - Ag1*(A11^-1*f1) - Ag2*(A22^-1*f2)
    for (int i=0; i<n; ++i) b[i] = fg[i] + grid1.v[n*i+n-1]*ih2 + grid2.v[n*i]*ih2;

    // Solve for ug using conjugate gradient
    conj_grad::solve();

    // Now use ug (stored in x) to solve on the subdomains
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            grid1.f[n*i+j] = f1[n*i+j];
            grid2.f[n*i+j] = f2[n*i+j];
        }
        grid1.f[n*i+n-1] += x[i]*ih2;
        grid2.f[n*i]     += x[i]*ih2;
    }
    grid1.solve(); // Compute u1 = A11^-1*(f1-A1g*ug)
    grid2.solve(); // Compute u2 = A22^-1*(f2-A2g*ug)
}

/** Perform matrix-vector multiplication with the Schur complement matrix S.
 * \param[in] in  : The input vector.
 * \param[in] out : The output vector, after multiplying by S. */
void schur::mul_A(double* in, double* out)
{
    // Use "in" as glue to solve on the subdomains
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            grid1.f[n*i+j] = 0;
            grid2.f[n*i+j] = 0;
        }
        grid1.f[n*i+n-1] = -in[i]*ih2;
        grid2.f[n*i]     = -in[i]*ih2;
    }

    grid1.solve(); // Compute A11^-1*A1g*in
    grid2.solve(); // Compute A22^-1*A2g*in

    // Compute out = S*in = Agg*in - Ag1*(A11^-1*A1g*in) - Ag2*(A22^-1*A2g*in)
    for (int i=0; i<n; ++i) {
        out[i] = 4*in[i] + grid1.v[n*i+n-1] + grid2.v[n*i];
        if (i>0)   out[i] -= in[i-1];
        if (i<n-1) out[i] -= in[i+1];
        out[i] *= ih2;
    }
}

/** Print the solution, padding with zeros for the Dirichlet boundary conditions. */
void schur::print_solution()
{
    for (int j=0; j<2*n+3; ++j) printf("0 ");
    puts("");
    for (int i=0; i<n; ++i) {
        printf("0 ");
        for (int j=0; j<n; ++j) printf("%g ", grid1.v[n*i+j]);
        printf("%g ", x[i]);
        for (int j=0; j<n; ++j) printf("%g ", grid2.v[n*i+j]);
        puts("0");
    }
    for (int j=0; j<2*n+3; ++j) printf("0 ");
}
