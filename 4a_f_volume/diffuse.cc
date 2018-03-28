#include <cstdlib>
#include <cstring>

#include "diffuse.hh"

/** Initializes the class for solving the diffusion equation on the periodic
 * unit interval, using a spatially dependent diffusion constant.
 * \param[in] m_ the number of gridpoints to use. */
diffuse::diffuse(int m_) : m(m_), dx(1./m), a(new double[m]), b(new double[m]),
        nu(new double[m]) {}

/** The class destructor frees the dynamically allocated memory. */
diffuse::~diffuse() {
    delete [] nu;
    delete [] b;
    delete [] a;
}

/** Performs one explicit timestep of the solution. The timestep itself it
 * already specified via the values in the nu array. */
template<int type>
void diffuse::step_forward() {
    for(int j=0;j<m;j++) {

        // Compute indices on left and right, taking into account periodicity
        int jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        // Perform update
        switch(type) {

            // Update type 0: conventional finite difference
            case 0: b[j]=((1-2*nu[j])*a[j]+nu[j]*(a[jl]+a[jr]));break;

            // Update type 1: finite volume method
            case 1: b[j]=a[j]+nu[j]*(a[jl]-a[j])+nu[jr]*(a[jr]-a[j]);
        }
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Solves the diffusion equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial snapshot).
 * \param[in] iters the number of iterations to step the solution forward by
 *                  between snapshots. */
void diffuse::solve(const char* filename,int snaps,int iters) {

    // Allocate memory to store solution snapshots. Integrate the system and
    // store snapshots after fixed time intervals
    double *z=new double[m*(snaps+1)];
    memcpy(z,a,m*sizeof(double));
    for(int i=1;i<=snaps;i++) {
        for(int k=0;k<iters;k++) step_forward();
        memcpy(z+i*m,a,m*sizeof(double));
    }

    // Open the output file to store the snapshots
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Can't open output file\n",stderr);
        exit(1);
    }

    // Print the snapshots, including periodic copies
    // at either end to get a full line over the interval
    // from 0 to 1
    print_line(-0.5*dx,z+(m-1),snaps);
    for(int j=0;j<m;j++) print_line((j+0.5)*dx,z+j,snaps);
    print_line(1+0.5*dx,z,snaps);

    // Delete snapshots and close file
    fclose(fp);
    delete [] z;
}

/** Initializes the solution to be a step function. */
void diffuse::init_step_function() {
    for(int i=0;i<m;i++) {
        x=dx*(i+0.5);
        a[i]=x>0.25&&x<0.75?1:0;
    }
}

/** Computes the integral of the solution over the domain.
 * \return The integral. */
double diffuse::integral() {
    double sum=0;
    for(double *ap=a;ap<a+m;ap++) sum+=*ap;
    return dx*sum;
}

/** Prints a line of stored snapshots to a file.
 * \param[in] fp a pointer to the file to write to.
 * \param[in] x the position in the domain corresponding to this line.
 * \param[in] zp a pointer to the first snapshot data point to print.
 * \param[in] snaps the number of snapshots (not including the starting
 *                  snapshot). */
void diffuse::print_line(FILE *fp,double x,double *zp,int snaps) {
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(" %g",zp[i*m]);
    fputc('\n',fp);
}

/** Initializes the nu array, which incorporates the timestep, grid spacing
 * and the spatially dependent diffusion constant.
 * \param[in] type the type of integration that will be performed. For type 0,
 *                 the nu values are stored at cell centers. For type 1, the nu
 *                 values are stored at cell edges.
 * \param[in] safe_fac a safety factor to multiply the maximally stable
 *                     timestep by. */
void diffuse::init_nu_array(int type,double safe_fac) {
    
    // Compute the timestep
    double pre=0.5*safe_fac/beta_max();
    dt=pre*dx*dx;

    // Compute the table of nu values
    double x=type==0?0.5*dx:0.;
    for(int i=0;i<m;i++,x+=dx) nu[i]=pre*beta(x);
}

// Explicit instantiation of the templated routine for the two different
// integration types
template void diffuse::step_forward<0>;
template void diffuse::step_forward<1>;
