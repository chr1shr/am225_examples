// DOP853 integrator re-written as a C++ class
//
// Author   : Chris H. Rycroft (Harvard SEAS / Lawrence Berkeley Laboratory)
// Email    : chr@alum.mit.edu
// Date     : May 6th 2015
//
// This source code file is based on the adaptive eighth-order Dormand-Prince
// integration code DOP853.F, developed by E. Hairer and G. Wanner. The
// original code is written is Fortran but has been refactored here into a C++
// class. For more information, see
//
// http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
//
// E. Hairer, S.P. Norsett AND G. Wanner, Solving ordinary differential
// equations I. Nonstiff problems, 2nd edition. Springer series in
// computational mathematics, Springer-Verlag (1993).

/** \file dop853.hh
 * \brief Header file for the dop853 class. */

#ifndef DOP853_HH
#define DOP853_HH

/** \brief Class for carrying out the adaptive eighth-order Dormand--Prince
 * integration scheme. */
class dop853 {
    public:
        /** The number of components in the differential equation system. */
        const int n;
        /** The maximum number of integration steps. */
        const long nmax;
        /** The number of timesteps between which to test for stiffness. */
        const long nstiff;
        /** A factor used in estimating the length of the next timestep. */
        const double fac1;
        /** A factor used in estimating the length of the next timestep. */
        const double fac2;
        /** A safety factor used in estimating the length of the next timestep.
         */
        const double safe;
        /** An exponent used in estimating the length of the next timestep. */
        const double beta;
        /** An estimate of the smallest number that can be added to 1 to obtain
         * a different number, used to determine when the timestep has become
         * too small. */
        const double uround;
        /** The absolute local error tolerance. */
        double atoli;
        /** The relative local error tolerance. */
        double rtoli;
        /** The total number of integration steps. */
        long nstep;
        /** The number of accepted integration steps. */
        long naccpt;
        /** The number of rejected integration steps. */
        long nrejct;
        /** The number of integration steps of fixed size. */
        long nfixed;
        /** The number of evaluations of the differential equation. */
        long nff;
        /** A counter for the number of non-stiff integration steps. */
        long nonsti;
        /** A counter for the number of stiff integration steps. */
        long iasti;
        /** A counter for the current snapshot number. */
        int csn;
        /** The current time. */
        double t;
        /** The previous time. */
        double told;
        /** The current time plus the integration step size. */
        double tph;
        /** The integration step size. */
        double h;
        /** An array holding the current state vector. */
        double *w;
        /** A temporary array for assembling a new state vector. */
        double *ww1;
        /** Storage for the Runge--Kutta intermediate steps. */
        double *k1,*k2,*k3,*k4,*k5,*k6,*k7,*k8,*k9,*k10;
        /** Storage for the dense output. */
        double *rc1,*rc2,*rc3,*rc4,*rc5,*rc6,*rc7,*rc8;
        dop853(int n_);
        ~dop853();
        int solve(double tstart,double tend,int snaps,double **ws);
        inline int solve(double tstart,double tend) {
            return solve(tstart,tend,0,NULL);
        }
        void fixed_step(double h_,bool last);
        bool detect_stiffness();
        void dense_output();
        void dense(double *ws,double ti);
        void init_counters();
        virtual void output() {};
        virtual void ff(double tt,double *in,double *out) = 0;
    private:
        void step12(double *p1);
        double error_estimation();
        double hinit(double hmax,double posneg);
        /** Returns the minimum of two numbers.
         * \param[in] (a,b) the two numbers.
         * \return The minimum */
        inline double min_d(double a,double b) {return a<b?a:b;}
        /** Returns the maximum of two numbers.
         * \param[in] (a,b) the two numbers.
         * \return The maximum */
        inline double max_d(double a,double b) {return a>b?a:b;}
};

#endif
