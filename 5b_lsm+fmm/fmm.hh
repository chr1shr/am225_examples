#ifndef FMM_HH
#define FMM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fields.hh"

const int 

/** A class to perform the fast marching method. */
class fmm {
    public:
        /** The number of grid cells in the horizontal direction. */
        const int m;
        /** The number of grid cells in the vertical direction. */
        const int n;
        /** The total number of grid cells. */
        const int mn;
        /** The memory step length, taking into account ghost point allocation. */
        const int ml;
        /** The periodicity in the x direction. */
        const bool x_prd;
        /** The periodicity in the y direction. */
        const bool y_prd;
        /** The lower bound in the x direction. */
        const double ax;
        /** The lower bound in the y direction. */
        const double ay;
        /** The upper bound in the x direction. */
        const double bx;
        /** The upper bound in the y direction. */
        const double by;
        /** The grid spacing in the x direction. */
        const double dx;
        /** The grid spacing in the y direction. */
        const double dy;
        /** The inverse grid spacing in the x direction. */
        const double xsp;
        /** The inverse grid spacing in the y direction. */
        const double ysp;
        /** An array containing the simulation fields. */
        phi_field* const phibase;
        /** A pointer to the (0,0) grid cell in the field array. */
        phi_field* const phim;
        fmm(const int m_,const int n_,const bool x_prd_,const bool y_prd_,
            const double ax_,const double bx_,const double ay_,const double by_);
        ~fmm();
    private:
	/** The number of elements on the heap. */
	int w;
	/** The current memory allocation for the heap. */
	int mem;
	/** The heap entries, which are pointers to the main grid. */
	phi_field** he;
        /** Temporary storage for used during the output routine. */
        float *buf;
};

#endif
