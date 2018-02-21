#ifndef SOL_RK4_HH
#define SOL_RK4_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class rk4 : public sol_base {
    public:
        rk4(int dof_);
        virtual ~rk4();
        virtual void step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
};

#endif
