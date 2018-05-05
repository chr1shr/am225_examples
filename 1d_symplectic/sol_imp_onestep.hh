#ifndef SOL_IMP_ONESTEP_HH
#define SOL_IMP_ONESTEP_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the one-step implicit method class,
 * which incorporates several different methods including the backward Euler
 * method, and the implicit midpoint method. */
class imp_onestep : public sol_base {
    public:
        /** The fraction of the timestep at which to compute the Runge--Kutta
         * step, which sets which method will be used. If alpha=0.5, this is
         * the implicit midpoint method. If alpha=1.0, this is the backward
         * Euler method. */
        const double alpha;
        imp_onestep(int dof_,double alpha_);
        virtual ~imp_onestep();
        virtual bool step(double dt);
    private:
        double *dq;
        double *k1;
        double *k1b;
};

#endif
