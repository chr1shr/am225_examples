#ifndef FLUID_2D_HH
#define FLUID_2D_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fields.hh"
#include "mgs_fem.hh"
#include "tgmg.hh"

/** A class to carry out a 2D incompressible fluid simulation. */
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
        int ntrace;
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
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        field* const fbase;
        /** A pointer to the (0,0) grid cell in the field array. */
        field* const fm;
        /** An array for the source terms used during the algebraic
         * multigrid solve. */
        double* const src;
        /** An array containing the tracer positions. */
        double* tm;
        /** The regular timestep to be used. */
        double dt_reg;
        /** The current simulation time. */
        double time;
        /** The current frame number. */
        int f_num;
        /** The files to save to file. */
        unsigned int fflags;
        fluid_2d(const int m_,const int n_,const bool x_prd_,const bool y_prd_,
                 const double ax_,const double bx_,const double ay_,const double by_,
                 const double visc_,const double rho_,unsigned int fflags_,
                 const char *filename_);
        ~fluid_2d();
        void solve(double duration,int frames);
        void step_forward(double dt);
        void init_fields();
        void write_files(int k);
        void initialize(int ntrace_,double dt_pad,double max_vel=-1);
        double advection_dt();
        void choose_dt(double dt_pad,double adv_dt,bool verbose=true);
        void init_tracers();
        void update_tracers(double dt);
        void output(const char *prefix,const int mode,const int sn,const bool ghost=false);
        void output_tracers(const char *prefix,const int sn);
        void save_header(double duration,int frames);
        /** Chooses a timestep size that is the largest value smaller than dt_reg,
        * such that a given interval length is a perfect multiple of this timestep.
        * \param[in] interval the interval length to consider.
        * \param[out] adt the timestep size.
        * \return The number of timesteps the fit into the interval. */
        inline int timestep_select(double interval, double &adt) {
            int l=static_cast<int>(interval/dt_reg)+1;
            adt=interval/l;
            return l;
        }
    private:
        /** An object containing the configuration of the first linear
         * system to be solved using the multigrid method. */
        mgs_fem ms_fem;
        void set_boundaries();
        void fem_source_term_conditions();
        double average_pressure();
        void copy_pressure();
        inline void vel_eno2(double &ud,double &vd,double hs,field &f0,field &f1,field &f2,field &f3);
        inline double eno2(double p0,double p1,double p2,double p3);
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
