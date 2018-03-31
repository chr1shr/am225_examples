#ifndef QUADRAT_HH
#define QUADRAT_HH

struct quadrat {
    /** The number of quadrature points. */
    const int n;
    /** The quadrature points, on the standard interval [-1,1]. */
    double* const x;
    /** The quadrature weights. */
    double* const w;
    quadrat(int n_);
    ~quadrat();
};

#endif
