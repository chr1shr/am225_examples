#ifndef LEVELSET_HH
#define LEVELSET_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "omp.h"

#include "fields.hh"

/** A class to carry out a simple level set simulation. */
class levelset {
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
        /** The square inverse grid spacing in the x direction. */
        const double xxsp;
        /** The square inverse grid spacing in the y direction. */
        const double yysp;
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        field* const fbase;
        /** A pointer to the (0,0) grid cell in the field array. */
        field* const fm;
        /** The regular timestep to be used. */
        double dt_reg;
        /** The current simulation time. */
        double time;
        /** The current frame number. */
        int f_num;
	/** Whether to oscillate the applied velocity field or not. */
	bool oscillate_vel;
        levelset(const int m_,const int n_,const bool x_prd_,const bool y_prd_,
                 const double ax_,const double bx_,const double ay_,const double by_,
                 const char *filename_);
        ~levelset();
        void solve(double duration,int frames);
        void step_forward(double dt);
        void init_fields(int type);
        inline void write_files(int k) {
            output(k,false);
        }
        void initialize(int type,bool oscillate_vel_,double dt_pad,double max_vel=-1);
        double advection_dt();
        void choose_dt(double dt_pad,double adv_dt,bool verbose=true);
        void output(const int sn,const bool ghost=false);
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
        void set_boundaries();
        inline void vel_eno2(double &phid,double hs,field &f0,field &f1,field &f2,field &f3);
        inline double eno2(double p0,double p1,double p2,double p3);
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
        /** Temporary storage for used during the output routine. */
        float *buf;
#ifdef _OPENMP
        inline double wtime() {return omp_get_wtime();}
#else
        inline double wtime() {return 0;}
#endif
};

#endif
