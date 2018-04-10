#include <cmath>

#include "cubic_1d_fe.hh"

int main() {

    // Construct the finite-element class and set the Neumann boundary
    // condition
    cubic_1d_fe cf(32);
    cf.g=exp(-1)*5*M_PI;

    // Initialize the source function to be a slope f(x)=1.5-x
    cf.init_mms();

    // Optional command to print the matrix in text form
    //cf.print_matrix();

    // Solve the finite-element problem using the conjugate gradient method
    cf.solve();
    cf.print();
}
