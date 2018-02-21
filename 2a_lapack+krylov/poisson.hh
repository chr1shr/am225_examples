#ifndef POISSON_HH
#define POISSON_HH

#include "conj_grad.hh"

class poisson : public conj_grad {
    public:
        /** The grid spacing. */
        const double dx;
        /** The inverse grid spacing squared. */
        const double xsp2;
        poisson(int dof_);
        void init();
        void print_solution();
        virtual void mul_A(double *in,double *out);
};

#endif
