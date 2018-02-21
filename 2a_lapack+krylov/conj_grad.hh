#ifndef CONJ_GRAD_HH
#define CONJ_GRAD_HH

#include "blas.h"

class conj_grad {
    public:
        const int dof;
        const int one;
        double *x;
        double *b;
        conj_grad(int dof_);
        ~conj_grad();
        void solve(bool verbose=false);
        void solve_pre(bool verbose=false);
        virtual void mul_A(double *in,double *out) = 0;
        virtual void M_inv(double *in,double *out);
    protected:
        inline void daxpy(double alpha,double *in,double *out) {
            daxpy_(&dof,&alpha,in,&one,out,&one);
        }
        inline void dscal(double alpha,double *x) {
            dscal_(&dof,&alpha,x,&one);
        }
        inline void copy(double *in,double *out) {
            dcopy_(&dof,in,&one,out,&one);
        }
        inline double dot(double *x,double *y) {
            return ddot_(&dof,x,&one,y,&one);
        }
    private:
        double *rk;
        double *pk;
        double *zk;
        double *yk;
};

#endif
