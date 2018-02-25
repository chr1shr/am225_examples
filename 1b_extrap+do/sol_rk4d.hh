#ifndef SOL_RK4D_HH
#define SOL_RK4D_HH

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class rk4d {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        rk4d(int dof_);
        virtual ~rk4d();
        void print(double t_,double *in);
        void dense_output(double theta,double dt);
        void solve_fixed(double t_end,int iters,bool output,int d_steps);
        void step(double dt);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
};

#endif
