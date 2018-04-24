#ifndef FIELDS_HH
#define FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
    /** The fixed horizontal velocity. */
	double u;
    /** The fixed vertical velocity. */
	double v;
    /** The level set function. */
    double phi;
    /** The change in the level set function during an update. */
    double cphi;
    inline void prd_bc(field &f) {
        phi=f.phi;
    }
    inline void neu_bc(field &f) {
        phi=f.phi;
    }
	inline void update() {
		phi+=cphi;
	}
    /** Computes the maximum allowable timestep based on the CFL restriction
     * from the velocity stored in this class.
     * \param[in] (xsp,ysp) the horizontal and vertical grid spacings.
     * \return The reciprocal of the maximum allowable timestep. */
    inline double cfl(double xsp,double ysp) {
        double uc=fabs(u)*xsp,vc=fabs(v)*ysp;
        return uc>vc?uc:vc;
    }
};

/** Data structure for the fast marching method example code. */
struct phi_field {
    /** The fast marching function. */
    double phi;
    /** An indicator variable. 0: empty, 1: on heap, 2: set, 3: boundary. */
    int c;
    /** A back pointer, indicating which position this gridpoint is on the
     * heap. */
    int bp;
};

#endif
