#include <cstdio>
#include <cmath>

int main() {
    double a[1024];

    // Since each entry of the array can be filled in separately, this
    // loop can be parallelized
#pragma omp parallel for
    for(int i=0;i<1024;i++) {
        a[i]=sqrt(double(i));
    }
}
