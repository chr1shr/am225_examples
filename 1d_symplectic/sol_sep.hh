#ifndef SOL_SEP_HH
#define SOL_SEP_HH

#include "sol_base.hh"

class sym_separable : public sol_base {
    public:
        /** Half the total degrees of freedom. */
        const int hdof;
        sym_separable(int dof_,int mode);
        virtual ~sym_separable();
        virtual bool step(double dt);
        virtual void ff(double t_,double *in,double *out) {}
        virtual void fp(double *in,double *out) = 0;
        virtual void fq(double *in,double *out) = 0;
    private:
        /** The number of steps in the method. */
        const int steps;
        /** The offset in the coefficient table for the method in use. */
        double* cstart;
        /** An array for computing the update to the momentum or coordinate
         * variables. */
        double* const del;
        static const int tab_steps[3];
        static const int tab_offset[3];
        static double coeffs[18];
};

#endif
