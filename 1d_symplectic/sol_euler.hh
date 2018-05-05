#ifndef SOL_EULER_HH
#define SOL_EULER_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the first-order explicit Euler method.
 */
class euler : public sol_base {
    public:
        euler(int dof_);
        virtual ~euler();
        virtual bool step(double dt);
    private:
        double *k1;
};

#endif
