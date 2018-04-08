#ifndef FLUID_2D_FIELDS_HH
#define FLUID_2D_FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /** The pressure. */
    double p;
    /** The change in horizontal velocity. */
    double cu;
    /** The change in vertical velocity. */
    double cv;
    inline void update() {
        u+=cu;
        v+=cv;
    }
};

#endif
