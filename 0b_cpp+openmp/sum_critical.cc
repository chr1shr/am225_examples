#include <cstdio>

int main() {
    unsigned int c=0;
#pragma omp parallel for
    for(unsigned int i=0;i<1024;i++) {
        int d=i*i;
#pragma omp critical
        {
            if(i%100==0) printf("Processing %d\n",i);
            c+=d;
        }
    }
    printf("Sum=%u\n",c);
}
