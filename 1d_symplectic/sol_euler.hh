#ifndef SOL_EULER_HH
#define SOL_EULER_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the first-order explicit Euler method.
 */
class euler : public sol_base {
    public:
        euler(int dof_);
        virtual ~euler();
        virtual void step(double dt);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        double *k1;
};

#endif
