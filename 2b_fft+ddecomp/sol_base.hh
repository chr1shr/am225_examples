#ifndef SOL_BASE_HH
#define SOL_BASE_HH

/** The base class for the ODE solvers, containing definitions and functions
 * that are common to all methods. */
class sol_base {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        sol_base(int dof_);
        virtual ~sol_base();
        void print();
        void solve_fixed(double t_end,int iters,bool output=false);
        virtual void step(double dt) = 0;
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
};

#endif
