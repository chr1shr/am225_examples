#include <cmath>
#include <cstdio>
#include <limits>

#include "quadrat.hh"

/** Initializes the Gaussian quadrature class. It computes the quadrature
 * points and weights by finding the roots of the associated Legendre
 * polynomial.
 * \param[in] n_ the number of quadrature points. */
quadrat::quadrat(int n_) : n(n_), x(new double[n]), w(new double[n]) {
    const double fac=M_PI/(0.5+n);
    double tmp,p1,p2,p3,pp,dz,z;
    int i,j,m=(n+1)>>1;

    // Due to symmetry, it is only necessary to consider the first half of the
    // points
    for(i=0;i<m;i++) {

        // Start with a guess for the root position
        z=cos((i+0.75)*fac);
        do {

            // Evaluate the Legendre polynomial at z using the recurrence relation
            p2=0.;
            p1=1.;
            for(j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1)*z*p2-(j-1)*p3)/static_cast<double>(j);
            }

            // Perform a Newton step
            tmp=n*(z*p1-p2);
            pp=tmp/(z*z-1.0);
            dz=-p1/pp;
            z+=dz;
        } while(fabs(dz)>=std::numeric_limits<double>::epsilon());

        // Store the quadrature point and weight. In addition, fill in the
        // (n-1-i)th point, using symmetry
        x[i]=-z;
        x[n-1-i]=-x[i];
        w[n-1-i]=w[i]=2.0/(-tmp*pp);
    }
}

/** The class destructor frees the dynamically allocated memory. */
quadrat::~quadrat() {
    delete [] w;
    delete [] x;
}
