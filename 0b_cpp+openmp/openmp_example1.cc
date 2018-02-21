#include <cstdio>

int main() {

#pragma omp parallel
    {
        // Since this is within a parallel block, each thread
        // will execute it
        puts("Hi");
    }
}
