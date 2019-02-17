#include <cstdio>
#include <cstdlib>

#include "lp_solve.hh"

// Generate uniform random numbers between -1 and 1
inline double rnd() {
    return -1+2./RAND_MAX*static_cast<double>(rand());
}

int main() {
    int i,j,k;

    // The size of the matrix
    const int n=6;

    // Number of lower and upper diagonals
    const int kl=2,ku=1;

    // Array for entries and source term
    const int ldab=2*kl+ku+1;
    double *A=new double[ldab*n],
           *b=new double[n];

    // Fill source term with random numbers
    for(j=0;j<n;j++) b[j]=rnd();

    // Fill matrix entries, putting ones on the diagonal, and small random
    // numbers in the off-diagonal locations
    for(i=0;i<n;i++) for(k=-kl;k<=ku;k++) {
        j=i-k;
        if(j<0||j>=n) continue;
        A[kl+ku+i*ldab-k]=k==0?1:0.1*rnd();
    }

    // Print the matrix
    for(j=0;j<n;j++) {
        printf(j==0?"A=[":"  [");
        for(i=0;i<n;i++) {
            k=i-j;
            printf("%7.2g ",k<-kl||k>ku?0:A[kl+ku+i*ldab-k]);
        }
        puts("]");
    }

    // Print the source term vector
    printf("\nb=[%7.3g ]\n",*b);
    for(j=1;j<n;j++) printf("  [%7.3g ]\n",b[j]);

    // Call routine to solve the matrix
    solve_banded_matrix(n,kl,ku,A,b);

    // Print the solution
    printf("\nx=[%7.3g ]\n",*b);
        for(j=1;j<n;j++) printf("  [%7.3g ]\n",b[j]);

    // Delete the dynamically allocated memory
    delete [] b;
    delete [] A;
}
