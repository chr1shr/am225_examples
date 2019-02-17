#include <cmath>
#include "schur.hh"

int main()
{
    // The number of units per subdomain in one dimension
    const int n = 100;
    auto f = [](double x, double y) { return y*y*std::exp(x-y); };
    schur poisson(n);
    poisson.solve(f);
    poisson.print_solution();
}
