#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "tgmg.hh"
#include "mgs_fem_vr.hh"

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif

int main(int argc,char **argv) {

    // Periodicity
    const bool x_prd=false,y_prd=false;

    // Grid dimensions
    const int m=129,n=129;

    // Physical grid size
    const double ax=-1,bx=1,ay=-1,by=1;

    // Grid spacings
    const double dx=(bx-ax)/(x_prd?m:m-1),
                 dy=(by-ay)/(y_prd?n:n-1);

    // Initialize the class in a test mode
    mgs_fem_varying_rho fvr(m,n,x_prd,y_prd,dx,dy);

    // Set up the 1/rho and source term fields
    int i,j,ij;
    double x,y;
    for(ij=j=0;j<n;j++) {
        y=ay+dy*j;
        for(i=0;i<m;i++,ij++) {
            x=ax+dx*i;

            // Set the 1/rho field to be a different value in a circle of
            // radius approximately 0.5.
            fvr.irho[ij]=x*x+y*y<0.2500123?3:1;

            // Set the source term. Note that since the BCs are either Neumann
            // (for non-periodic), or are periodic, it is necessary that the
            // source term intergrate to zero.
            fvr.src[ij]=(x-0.7)*(x-0.7)+(y-0.7)*(y-0.7)<0.04000789?dx*dy:
                       ((x+0.7)*(x+0.7)+(y+0.7)*(y+0.7)<0.04000789?-dx*dy:0);
        }
    }

    // Setup the linear systems on the coarser levels of the multigrid
    // hierarchy. IMPORTANT NOTE: when the irho field changes, this function
    // must be called prior to performing V-cycles, since the coarse linear
    // systems will change
    fvr.mg.verbose=3;
    double t0=wtime(),t1;
    fvr.setup();

    // Perform V-cycles and print out timing information
    t1=wtime();
    fvr.solve_v_cycle();
    printf("# Setup time : %g ms\n"
           "# Solve time : %g ms\n",1e3*(t1-t0),1e3*(wtime()-t1));

    // Save the source term and solution in a format that can be read by Gnuplot.
    fvr.mg.output_b("b.0",ax,dx,ay,dy);
    fvr.mg.output_z("z.0",ax,dx,ay,dy);
}
