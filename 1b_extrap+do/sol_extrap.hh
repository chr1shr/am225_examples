#ifndef SOL_EXTRAP_HH
#define SOL_EXTRAP_HH

/** Class for solving an ODE IVP using the extrapolation procedure. */
class extrap {
    public:
        /** Whether to use the Gragg method for basic integration, which
         * doubles the order of accuracy. */
        const bool gragg;
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The maximum allowed value of the j index in the extrapolation
         * method. */
        int max_j;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        extrap(bool gragg_,int dof_,int max_j_);
        virtual ~extrap();
        void init_romberg();
        void init_bulirsch();
        void init_harmonic();
        void print(double t_,double *in);
        void solve_fixed(double t_end,int iters,bool output,int j,int k);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        void compute_basic(int j,double t);
        double* dq;
        double* k1;
        double* aq;
        int* n_seq;
};

#endif
