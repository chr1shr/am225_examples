#include <cstdio>
#include <cstring>

int main() {

    // Grid size
    const int m=1024;
    double a[m],b[m],x;

    // Number of snapshots to output, and iterations between snapshot
    const int snaps=40,iters=20000;
    
    // Snapshot storage
    double z[m*(snaps+1)];

    // PDE-related constants
    const double alpha=0.1,
                 dx=1.0/m,
                 dt=0.1*dx*dx/alpha,
                 nu=alpha*dt/(dx*dx);

    // Initial condition
    for(int i=0;i<m;i++) {
        x=dx*i;
        a[i]=x>0.25&&x<0.75?1:0;
    }

    // Use system function to copy initial condition into snapshot storage
    memcpy(z,a,m*sizeof(double));

    // Integrate the PDE
    for(int i=1;i<=snaps;i++) {
        for(int k=0;k<iters;k++) {

            // This loop can be parallelized using OpenMP, by uncommenting the
            // line below
//#pragma omp parallel for
            for(int j=0;j<m;j++) {
                int jl,jr;
            
                // Compute indices on left and right, taking into account periodicity
                jl=j==0?m-1+j:j-1;
                jr=j==m-1?1-m+j:j+1;
                
                // Perform update
                b[j]=((1-2*nu)*a[j]+nu*(a[jl]+a[jr]));
            }
            memcpy(a,b,m*sizeof(double));
        }
        memcpy(z+i*m,a,m*sizeof(double));
    }

    // Output the results
    for(int j=0;j<m;j++) {
        printf("%g",j*dx);
        for(int i=0;i<=snaps;i++) printf(" %g",z[j+i*m]);
        putc('\n',stdout);
    }
}
