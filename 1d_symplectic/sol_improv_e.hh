#ifndef SOL_IMPROV_E_HH
#define SOL_IMPROV_E_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the improved Euler second-order method.
 */
class improv_e : public sol_base {
    public:
        improv_e();
        improv_e(int dof_);
        virtual ~improv_e();
        virtual bool step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
};

#endif
