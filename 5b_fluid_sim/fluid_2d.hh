#ifndef FLUID_2D_HH
#define FLUID_2D_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fields.hh"
#include "multisetup.hh"
#include "tgmg.hh"
#include "visco_impl.hh"

/** \brief A class to carry out a 2D elasticity simulation. */
class fluid_2d {
    public:
        /** The number of grid cells in the horizontal direction. */
        const int m;
        /** The number of grid cells in the vertical direction. */
        const int n;
        /** The total number of grid cells. */
        const int mn;
        /** The number of cell-cornered points in the pressure
         * finite-element problem in the horizontal direction. */
        const int m_fem;
        /** The number of cell-cornered points in the pressure
         * finite-element problem in the vertical direction. */
        const int n_fem;
        /** The memory step length, taking into account ghost point allocation. */
        const int ml;
        /** The number of tracers. */
        const int ntrace;
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
        /** The square inverse grid spacing in the x direction. */
        const double xxsp;
        /** The square inverse grid spacing in the y direction. */
        const double yysp;
        /** The viscosity. */
        const double visc;
        /** The density. */
        const double rho;
        /** The inverse density. */
        const double rhoinv;
        /** A multiplier to apply to the default timestep size. */
        const double tmult;
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        field* const fbase;
        /** A pointer to the (0,0) grid cell in the field array. */
        field* const fm;
        /** An array containing the tracer positions. */
        double* const tm;
        /** An array for the source terms used during the algebraic
         * multigrid solve.*/
        double* const src;
        /** An array for the finite-element solution. */
        double* const sfem;
        /** The current simulation time. */
        double time;
        fluid_2d(const int m_,const int n_,const bool x_prd,const bool y_prd,
             const double ax_,const double bx_,const double ay_,const double by_,
             const double re_,const double rho_,const double tmult_,
             const int ntrace_,const bool implicit_visc,const char *filename_);
        ~fluid_2d();
        void solve(double t_start,double t_end,int frames);
        void solve(double dt, int frames, bool pres, bool visc, bool force);
        void step_forward(double dt);
        void init_fields();
        void write_files(int k);
        void write_files(const char *fn, int k);
        void init_tracers();
        void update_tracers(double dt);
        void output(const char *prefix,const int mode,const int sn,const bool ghost=false);
        void output_tracers(const char *prefix,const int sn);
    private:
        /** An object containing the configuration of the first linear
         * system to be solved using the multigrid method. */
        msu_fem m_fem;
        /** The multigrid solver for the first multigrid problem. */
        tgmg<msu_fem,double,double> m_fem;
        void set_boundaries();
        void fem_source_term_conditions();
        double average_pressure();
        double mymax(double *G, int nib);
        double mymin(double *G, int nib);
        void copy_pressure();
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
        inline void remap_tracer(double &xx,double &yy);
        /** Temporary storage for used during the output routine. */
        float *buf;
#ifdef _OPENMP
        inline double wtime() {return omp_get_wtime();}
#else
        inline double wtime() {return 0;}
#endif
};

#endif
