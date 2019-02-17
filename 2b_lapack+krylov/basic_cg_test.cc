#include "poisson.hh"

int main() {
    poisson po(11);
    po.init();
    po.solve(true);
    po.print_solution();
}
