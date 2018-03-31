#include <cstdio>
#include <cmath>

#include "quadrat.hh"

// Two-dimensional test function
double f(double x,double y) {
    return cos(M_PI*(x+0.25*y))*exp(-x);
}

// Exact integral of test function on the square [-1,1]^2, calculated using
// Mathematica
const double exact=4*sqrt(2)*(1-exp(2))/(M_PI*exp(1)*(1+M_PI*M_PI));

int main() {
    int i,j,n;
    double I,Ir;

    // Loop over a range of numbers of quadrature points
    for(n=1;n<=12;n++) {

        // Set up the quadrature points and weights
        quadrat q(n);

        // Perform the 2D sum of function evaluations, each multiplied by the
        // corresponding weight
        I=0.;
        for(j=0;j<n;j++) {
            Ir=0.;
            for(i=0;i<n;i++) Ir+=q.w[i]*f(q.x[i],q.x[j]);
            I+=Ir*q.w[j];
        }

        // Print the integral value and the absolute error from the exact
        // solution
        printf("%d %g %g\n",n,I,I-exact);
    }
}
