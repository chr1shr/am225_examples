#ifndef FIELDS_HH
#define FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
	double u;
	double v;
    double phi;
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

#endif
