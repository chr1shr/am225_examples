#include "schwarz.hh"

#include <cstdio>
#include "omp.h"

int main() {
    double t0,t1;

    schwarz sz(64,32,32);

    sz.init();
    t0=omp_get_wtime();
    sz.solve_additive();
    t0=omp_get_wtime()-t0;
    putchar('\n');

    sz.clear_solution();
    t1=omp_get_wtime();
    sz.solve_multiplicative();
    t1=omp_get_wtime()-t1;
    putchar('\n');

    sz.output_solution("add_sch.sol");
    sz.output_source("add_sch.src");

    printf("# Additive: %g s\n# Multiplicative %g s\n",t0,t1);
}
