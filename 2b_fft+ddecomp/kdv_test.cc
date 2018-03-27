#include <cstdio>
#include <cmath>

#include "kdv.hh"

const int n=128;

int main() {
    int j,l;
    kdv k(n);
    double h=2*M_PI/n;

    k.init();
    for(int i=0;i<n;i++) printf("%g %g\n",h*i,k.q[i]);
    puts("\n");

    for(j=0;j<100;j++) {
        for(l=0;l<50;l++) k.step(0.001);
        for(int i=0;i<n;i++) printf("%g %g\n",h*i,k.q[i]);
        puts("\n");
    }
}
