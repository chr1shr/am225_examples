#ifndef SOL_HEUN3_HH
#define SOL_HEUN3_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the third-order Heun method. */
class heun3 : public sol_base {
    public:
        heun3();
        heun3(int dof_);
        virtual ~heun3();
        virtual void step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
};

#endif
