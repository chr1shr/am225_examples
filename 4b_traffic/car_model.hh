#ifndef CAR_MODEL_HH
#define CAR_MODEL_HH

#include "sol_rk4.hh"

/** This class has functions to specify the discrete car model. */
class car_model : public rk4 {
    public:
        /** The length of the periodic domain. */
        const double L;
        car_model(int dof_) : rk4(dof_), L(1.5*dof_) {};
        virtual void ff(double t_,double *in,double *out);
        virtual void init();
};

#endif
