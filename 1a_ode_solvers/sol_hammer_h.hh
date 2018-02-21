#ifndef SOL_HAMMER_H_HH
#define SOL_HAMMER_H_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class hammer_h : public sol_base {
    public:
        hammer_h(int dof_);
        virtual ~hammer_h();
        virtual void step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k1b;
        double *k2b;
};

#endif
