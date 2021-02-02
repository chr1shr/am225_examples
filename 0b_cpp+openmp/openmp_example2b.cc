#include <cstdio>

// A major feature of OpenMP is that if the code is compiled without OpenMP
// support, then the #pragma omp commands are ignored and the code reduces to a
// regular C++ program. This is very helpful and allows compatibility with
// systems where OpenMP is unavailable.
//
// An exception to this is if the "omp.h" header is included to access the
// OpenMP-specific functions. These functions won't be available without
// OpenMP. However, a program can still be made compatible with non-OpenMP
// compilation. The preprocessor variable _OPENMP is only defined when
// OpenMP is available. This variant of openmp_example2.cc shows a way
// to wrap the OpenMP functions, and provide alternatives in the case
// when OpenMP isn't available.

#ifdef _OPENMP

// OpenMP is available. Define thread_num and max_threads to call the
// corresponding OpenMP routines.
#include "omp.h"
inline int thread_num() {
    return omp_get_thread_num();
}
inline int max_threads() {
    return omp_get_max_threads();
}
#else

// OpenMP is not available. Define thread_num and max_threads to just return
// some simple values consistent with running on a single thread.
inline int thread_num() {return 0;}
inline int max_threads() {return 1;}
#endif

int main() {

#pragma omp parallel
    {
        // Variables declared within a parallel block are local to it
        int i=thread_num(),j=max_threads();

        printf("Hello from thread %d of %d\n",i,j);
    }
}
