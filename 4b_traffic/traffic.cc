#include <cstdlib>
#include <cstring>

#include "traffic.hh"

/** Initializes the class for solving the traffic equation on the periodic unit
 * interval.
 * \param[in] m_ the number of gridpoints to use. */
traffic::traffic(int m_,double A_) : m(m_), dx(1./m),
    a(new double[m]), b(new double[m]), F(new double[m]) {}

/** The class destructor frees the dynamically allocated memory. */
traffic::~traffic() {
    delete [] F;
    delete [] b;
    delete [] a;
}

/** Performs one explicit timestep of the nonlinear Godunov method.
 * \param[in] dt the timestep to use. */
void traffic::godunov(double dt) {
    double f=dt/dx,s;

    for(int j=0;j<m;j++) {
        int jl=j==0?m-1+j:j-1;

        if(a[j]<0.5&&a[jl]>=0.5) {
            F[j]=0.25;
        } else {
            s=(flux(a[j])-flux(a[jl]))/(a[j]-a[jl]);
            F[j]=s<0?flux(a[j]):flux(a[jl]);
        }
    }

    for(int j=0;j<m;j++) {

        // Compute index to the right and perform update
        int jr=j==m-1?1-m+j:j+1;
        b[j]=a[j]-f*(F[jr]-F[j]);
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Solves the transport equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial
 *                  snapshot).
 * \param[in] duration the number of iterations to step the solution forward by
 *                     between snapshots.
 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction. */
void traffic::solve(const char* filename,int snaps,double duration,double safe_fac) {

    // Compute the timestep and number of iterations
    double interval=duration/snaps,dt=dx*safe_fac;
    int iters=static_cast<int>(interval/dt)+1;
    dt=interval/iters;

    // Allocate memory to store solution snapshots. Integrate the system and
    // store snapshots after fixed time intervals
    double *z=new double[m*(snaps+1)];
    memcpy(z,a,m*sizeof(double));
    for(int i=1;i<=snaps;i++) {

        // Perform the explict timesteps
        for(int k=0;k<iters;k++) godunov(dt);

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
void traffic::init_step_function() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=x>0.25&&x<0.75?0.8:0;
    }
}

/** Initializes the solution to be the exponential of a sine wave. */
void traffic::init_exp_sine() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=exp(sin(2*M_PI*x));
    }
}

/** Computes the integral of the solution over the domain.
 * \return The integral. */
double traffic::integral() {
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
void traffic::print_line(FILE *fp,double x,double *zp,int snaps) {
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",zp[i*m]);
    fputc('\n',fp);
}
