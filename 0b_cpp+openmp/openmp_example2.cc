#include <cstdio>

// OpenMP header file with specific thread-related functions
#include "omp.h"

int main() {

#pragma omp parallel
    {
        // Variables declared within a parallel block are local to it
        int i=omp_get_thread_num(),j=omp_get_max_threads();

        printf("Hello from thread %d of %d\n",i,j);
    }
}
