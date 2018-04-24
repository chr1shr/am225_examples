#ifndef FMM_HH
#define FMM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fields.hh"

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
        /** The square inverse grid spacing in the x direction. */
        const double xxsp;
        /** The square inverse grid spacing in the y direction. */
        const double yysp;
        /** An array containing the simulation fields. */
        phi_field* const phibase;
        /** A pointer to the (0,0) grid cell in the field array. */
        phi_field* const phim;
        fmm(const int m_,const int n_,const double ax_,const double bx_,
            const double ay_,const double by_);
        ~fmm();
        void init_fields();
        void mark_box(int il,int iu,int jl,int ju);
        void init_heap();
        void fast_march();
        void add_neighbors(phi_field *phip);
        void add_to_heap(phi_field *phip);
        double calc_phi(phi_field *phip);
        double phi_full(double phiv,double phih);
        void setup_indicator_field();
        void trickle(phi_field *phip);
        void update(phi_field *phip);
        void reduce_heap();
        void set_boundary_phi();
        void integrate_path(double x,double y);
        double bilinear(double x,double y,double &gx,double &gy);
        void output(const char *filename,int mode);
    private:
        inline double min(double a,double b) {
            return a<b?a:b;
        }
        void add_heap_memory();
        inline bool phi_look(phi_field *phip,int d,double &phid);
        /** The heap counter, equal to the number of elements minus one. */
        int w;
        /** The current memory allocation for the heap. */
        int mem;
        /** The heap entries, which are pointers to the main grid. */
        phi_field** he;
        /** Temporary storage for used during the output routine. */
        float *buf;
};

#endif
