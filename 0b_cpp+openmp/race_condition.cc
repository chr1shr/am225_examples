#include <cstdio>
    
int main() {
    unsigned int c=0;
#pragma omp parallel for
    for(unsigned int i=0;i<1024;i++) {
        c=i*i+c;
    }
    printf("Sum=%u\n",c);
}
