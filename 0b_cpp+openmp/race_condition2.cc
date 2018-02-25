#include <cstdio>
    
int main() {
    int c[4096],d;

    // Fill table with square numbers
#pragma omp parallel for
    for(int i=0;i<4096;i++) {
        d=i*i;
        c[i]=d;
    }

    // Print out discrepancies
    for(int i=0;i<4096;i++)
    if(c[i]!=i*i) printf("%d %d\n",i,c[i]);
}
