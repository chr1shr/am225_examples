#include <cstdio>
#include <cstring>
#include <cmath>

#include "omp.h"

double beta(double x) {
	return 0.12+0.08*sin(x);
}

int main() {

    // Grid size
    const int m=512;
    double a[m],b[m],x;

    // Number of snapshots to output, and iterations between snapshot
    const int snaps=40,iters=20000;

    // Snapshot storage
    double z[m*(snaps+1)];

    // PDE-related constants
    const double alpha=0.2,
                 dx=1.0/m,
                 dt=0.1*dx*dx/alpha;
    double nu;

    // Initial condition
    for(int i=0;i<m;i++) {
        x=dx*i;
        a[i]=x>0.25&&x<0.75?1:0;
    }

    // Use system function to copy initial condition into snapshot storage
    memcpy(z,a,m*sizeof(double));

    // Integrate the PDE
    double t0=omp_get_wtime();
    for(int i=1;i<=snaps;i++) {
        for(int k=0;k<iters;k++) {
            for(int j=0;j<m;j++) {
                int jl,jr;

                // Compute indices on left and right, taking into account periodicity
                jl=j==0?m-1+j:j-1;
                jr=j==m-1?1-m+j:j+1;

                // Perform update
		nu=beta(j*dx);
                //b[j]=((1-2*nu)*a[j]+nu*(a[jl]+a[jr]));

		// Perform update 2
		b[j]=a[j]+beta((j-0.5)*dx)*(a[jl]-a[j])+beta((j+0.5)*dx)*(a[jr]-a[j]);
            }
            memcpy(a,b,m*sizeof(double));
	    
        }
//	double bs=0;
//	for(int j=0;j<m;j++) bs+=b[j];
//	printf("%d %g\n",i,bs);
        memcpy(z+i*m,a,m*sizeof(double));
    }
//    printf("# Wall clock time: %g s\n",omp_get_wtime()-t0);

    // Output the results
    for(int j=0;j<m;j++) {
        printf("%g",j*dx);
        for(int i=0;i<=snaps;i++) printf(" %g",z[j+i*m]);
        putc('\n',stdout);
    }
}
