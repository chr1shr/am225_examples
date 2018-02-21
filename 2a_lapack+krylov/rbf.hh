#ifndef RBF_HH
#define RBF_HH

#include <cstdlib>

#include "conj_grad.hh"

class rbf : public conj_grad {
    public:
        /** Number of specified points. */
        const int n;
        /** The type of radial function to use. */
        const int type;
        /** The square length scale of the radial function. */
        double lsq;
        /** The inverse square length scale of the radial function. */
        double ilsq;
        /** x positions of points. */
        double* const px;
        /** y positions of points. */
        double* const py;
        /** Function values at points. */
        double* const pf;
        /** Weights of the points in the RBF interpolant. */
        double* const rs;
        rbf(int n_,int type_);
        ~rbf();
        void init_random();
        void eigenvalues();
        void solve_weights_lapack();
        void solve_weights_conj_grad(int bls=0,bool verbose=false);
        void output_points(const char* filename);
        void output_interpolant(const char* filename,int q,double L);
        virtual void mul_A(double *in,double *out);
        virtual void M_inv(double *in,double *out);
        void make_table();
        inline void set_length_scale(double lscale) {
            lsq=lscale*lscale;ilsq=1./lsq;
        }
        inline void assemble_matrix() {
            if(A==NULL) fill_matrix_entries(A=new double[n*n],0,n);
        }
    protected:
        double urand();
    private:
        void preconditioning_table(int bls_);
        void fill_matrix_entries(double *Ap,int k,int b);
        double phi(double rsq);
        /** The block size in the block Jacobi preconditioner. */
        int bls;
        /** The last block size in the block Jacobi preconditioner. */
        int lbls;
        /** A pointer to the dense matrix entries. */
        double *A;
        /** A pointer to the block Jacobi preconditioning information. */
        double *Apre;
        /** A pointer to the pivoting information in the block Jacobi
         * preconditioner. */
        int *ipiv;
        /** The number of entries in each row of the sparse matrix. */
        int* An;
        /** The indices of non-zero entries in the sparse matrix. */
        int* phi_index;
        /** The values of non-zero entries in the sparse matrix. */
        double* phi_val;
};

#endif
