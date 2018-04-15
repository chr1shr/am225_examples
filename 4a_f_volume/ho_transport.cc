#include <cstdlib>
#include <cstring>

#include "ho_transport.hh"

/** Initializes the class for solving the transport equation on the periodic
 * unit interval, using a variety of high-order and high-resolution methods.
 * \param[in] m_ the number of gridpoints to use.
 * \param[in] A_ the advection velocity. */
ho_transport::ho_transport(int m_,double A_) : m(m_), dx(1./m), A(A_),
    a(new double[m]), b(new double[m]), sdx(new double[m]) {}

/** The class destructor frees the dynamically allocated memory. */
ho_transport::~ho_transport() {
    delete [] sdx;
    delete [] b;
    delete [] a;
}

/** Performs one step of the first-order Godunov method using a piecewise
 * constant reconstruction. This method is equivalent to a one-sided
 * first-order finite difference update, and is provided for comparison
 * purposes.
 * \param[in] dt the timestep to use. */
void ho_transport::godunov(double dt) {
    double f=A/dx*dt;

    // Perform one-sided update, taking into account periodicity
    for(int j=0;j<m;j++)
        b[j]=a[j]+f*(a[j==0?m-1+j:j-1]-a[j]);

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Performs one explicit timestep of the Lax-Wendroff method.
 * \param[in] dt the timestep to use. */
void ho_transport::lax_wendroff(double dt) {
    double f=A/dx*dt,
           sl=0.5*f*(1+f),sc=1-f*f,sr=0.5*f*(-1+f);

    for(int j=0;j<m;j++) {

        // Compute indices on left and right, taking into account periodicity
        int jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        // Perform update
        b[j]=sl*a[jl]+sc*a[j]+sr*a[jr];
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Performs one explicit timestep of the Beam-Warming method.
 * \param[in] dt the timestep to use. */
void ho_transport::beam_warming(double dt) {
    double f=A/dx*dt,
           sll=-0.5*f+0.5*f*f,sl=f*(2-f),sc=1-1.5*f+0.5*f*f;

    for(int j=0;j<m;j++) {

        // Compute two indices to the left, taking into account periodicity
        int jl=j==0?m-1+j:j-1,
            jll=j<=1?m-2+j:j-2;

        // Perform update
        b[j]=sll*a[jll]+sl*a[jl]+sc*a[j];
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Performs one explicit timestep using a slope-limited method.
 * \param[in] dt the timestep to use.
 * \param[in] type the integration type to use. 2:minmod, 3:superbee. */
void ho_transport::slope_limiter(double dt,int type) {

    // Set up slopes using either the minmod or superbee slope limiting
    // procedure
    type==3?sl_setup<3>():sl_setup<4>();

    double f=A/dx*dt;
    for(int j=0;j<m;j++) {

        // Compute index on left, taking into account periodicity
        int jl=j==0?m-1+j:j-1;

        // Perform update
        b[j]=(1-f)*a[j]+f*a[jl]-0.5*f*(1-f)*(sdx[j]-sdx[jl]);
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Fills the (sigma*dx) slope array using a limiting procedure. If type=3, use
 * the minmod method. If type=4, use the superbee method. */
template<int type>
void ho_transport::sl_setup() {
    for(int j=0;j<m;j++) {

        // Compute indices on left and right, taking into account periodicity
        int jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        // Apply the slope limiter
        sdx[j]=type==3?minmod(a[j]-a[jl],a[jr]-a[j])
                      :maxmod(minmod(a[jr]-a[j],2*(a[j]-a[jl])),
                              minmod(2*(a[jr]-a[j]),a[j]-a[jl]));
    }
}

/** Solves the transport equation using the second-order essentially
 * non-oscillatory (ENO) method. Note this is not a finite-volume method (and
 * is not conservative) but it still does a good job at dealing with sharp
 * features in a solution.
 * \param[in] dt the timestep to use. */
void ho_transport::step_eno2(double dt) {
    double f=0.5*A/dx*dt;
    for(int j=0;j<m;j++) {

        // Compute required indices, taking into account periodicity
        int jl=j==0?m-1+j:j-1,jll=j<=1?m-2+j:j-2,
            jr=j==m-1?1-m+j:j+1;

        // Perform update
        b[j]=a[j]-f*eno2(a[jr],a[j],a[jl],a[jll]);
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Calculates the ENO derivative using a sequence of values at four
 * gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double ho_transport::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Solves the transport equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial
 *                  snapshot).
 * \param[in] duration the number of iterations to step the solution forward by
 *                     between snapshots.
 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction.
 * \param[in] type the integration type to use. 0: Lax-Wendroff,
 *                 1: Beam-Warming, 2: minmod, 3: superbee. */
void ho_transport::solve(const char* filename,int snaps,double duration,double safe_fac,int type) {

    // Compute the timestep and number of iterations
    double interval=duration/snaps,dt=dx/A*safe_fac;
    int iters=static_cast<int>(interval/dt)+1;
    dt=interval/iters;

    // Allocate memory to store solution snapshots. Integrate the system and
    // store snapshots after fixed time intervals
    double *z=new double[m*(snaps+1)];
    memcpy(z,a,m*sizeof(double));
    for(int i=1;i<=snaps;i++) {

        // Perform the explict timesteps
        switch(type) {
            case 0: for(int k=0;k<iters;k++) godunov(dt);break;
            case 1: for(int k=0;k<iters;k++) lax_wendroff(dt);break;
            case 2: for(int k=0;k<iters;k++) beam_warming(dt);break;
            case 3:
            case 4: for(int k=0;k<iters;k++) slope_limiter(dt,type);break;
            case 5: for(int k=0;k<iters;k++) step_eno2(dt);
        }

        // Store the snapshot
        memcpy(z+i*m,a,m*sizeof(double));
    }

    // Open the output file to store the snapshots
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Can't open output file\n",stderr);
        exit(1);
    }

    // Print the snapshots, including periodic copies at either end to get a
    // full line over the interval from 0 to 1
    print_line(fp,-0.5*dx,z+(m-1),snaps);
    for(int j=0;j<m;j++) print_line(fp,(j+0.5)*dx,z+j,snaps);
    print_line(fp,1+0.5*dx,z,snaps);

    // Delete snapshots and close file
    fclose(fp);
    delete [] z;
}

/** Initializes the solution to be a step function. */
void ho_transport::init_step_function() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=x>0.1&&x<0.6?1:0;
    }
}

/** Initializes the solution to be the exponential of a sine wave. */
void ho_transport::init_exp_sine() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=exp(sin(2*M_PI*x));
    }
}

/** Computes the integral of the solution over the domain.
 * \return The integral. */
double ho_transport::integral() {
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
void ho_transport::print_line(FILE *fp,double x,double *zp,int snaps) {
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",zp[i*m]);
    fputc('\n',fp);
}

// Explicit instantiation
template void ho_transport::sl_setup<3>();
template void ho_transport::sl_setup<4>();
