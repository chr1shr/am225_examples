#ifndef SOL_RALSTON_HH
#define SOL_RALSTON_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using Ralston's second-order method. */
class ralston : public sol_base {
    public:
        ralston();
        ralston(int dof_);
        virtual ~ralston();
        virtual void step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
};

#endif
