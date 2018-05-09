#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "v_visco.hh"

/** Sets up the vanishing viscosity class, initializing constants and
 * allocating memory.
 * \param[in] m_ the number of gridpoints.
 * \param[in] epsilon_ the diffusion constant. */
v_visco::v_visco(int m_,double epsilon_) : m(m_), dx(1./m), epsilon(epsilon_),
    a(new double[m]), b(new double[m]), z(NULL), snaps(0) {}

/** The class destructor frees the dynamically allocated memory. */
v_visco::~v_visco() {
    if(z!=NULL) delete [] z;
    delete [] b;
    delete [] a;
}

/** Sets up the initial condition. */
void v_visco::initialize() {
    double fac=2*M_PI*dx;
    for(int i=0;i<m;i++) a[i]=exp(-2-2*cos(fac*i));//+0.4*sin(fac*i);
    //for(int i=0;i<m;i++) a[i]=0.5+0.4*sin(fac*i);
}

/** Integrates the transport equation using the vanishing viscosity approach
 * over a given duration.
 * \param[in] duration the duration to integrate over.
 * \param[in] sfac a safety factor to multiply the maximum stable timestep by.
 * \param[in] save whether to save snapshots of the simulation at periodic
 *                 intervals.
 * \param[in] snaps_ the number of snapshots to save (only used if 'save' is
 *                   set to true). */
void v_visco::integrate(double duration,double sfac,bool save,int snaps_) {

    // For cases with no output, overwrite the snapshot number to 1 to get
    // correct timestep calculation
    if(!save) snaps_=1;

    // Compute initial timestep, and time between snapshots
    double dt_base=sfac*dx*dx/epsilon,
           snap_dur=duration/snaps_,dt_adv=dx;
    if(dt_adv<dt_base) dt_base=dt_adv;

    // Adjust timestep down, so that an exact number of timesteps will fit in
    // the snapshot duration
    int iters=int(snap_dur/dt_base)+1;
    double dt=snap_dur/iters;

    // If snapshot saving is requested, then set up memory and save the initial
    // condition
    if(save) {
        if(z!=NULL) delete [] z;
        snaps=snaps_;
        z=new double[m*(snaps+1)];
        memcpy(z,a,m*sizeof(double));
    }

    // Integrate the PDE
    for(int i=1;i<=snaps_;i++) {
        for(int k=0;k<iters;k++) step(dt);
        if(save) memcpy(z+i*m,a,m*sizeof(double));
    }
}

/** Performs an explicit finite-difference timestep update of the traffic equation
 * using a small added viscosity.
 * \param[in] dt the timestep to take. */
void v_visco::step(double dt) {
    double f=0.5*dt/dx,nu=epsilon*dt/(dx*dx);
    for(int j=0;j<m;j++) {

        // Compute indices on left and right, taking into account periodicity
        int jll=j<=1?m-2+j:j-2,
            jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        // Calculate the finite-difference update at this gridpoint
        b[j]=(1-2*nu)*a[j]+nu*(a[jl]+a[jr])-f*(1-2*a[j])*eno2(a[jr],a[j],a[jl],a[jll]);
    }

    // Copy the new values over the old
    memcpy(a,b,m*sizeof(double));
}

/** Calculates the ENO derivative using a sequence of values at four
 * gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double v_visco::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Outputs the previously stored snapshots.
 * \param[in] filename the name of the file to write to. */
void v_visco::output(const char* filename) {

    // Open the output file in write mode "w", and give an error if there is a
    // failure (signified by a NULL pointer on return).
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening output file\n",stderr);
        exit(1);
    }

    // Check that the snapshots are actually allocated
    if(z==NULL) {
        fputs("Snapshots were not previously allocated\n",stderr);
        exit(1);
    }

    // Output the results
    for(int j=0;j<m;j++) {
        fprintf(fp,"%g",j*dx);
        for(int i=0;i<=snaps;i++) fprintf(fp," %g",z[j+i*m]);
        putc('\n',fp);
    }

    // Output the final line
    fputs("1",fp);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",z[i*m]);
    putc('\n',fp);
    fclose(fp);
}
