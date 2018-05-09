#ifndef V_VISCO_HH
#define V_VISCO_HH

class v_visco {
    public:
        /** The number of grid points. */
        const int m;
        /** The grid spacing. */
        const double dx;
        /** The viscosity. */
        const double epsilon;
        /** The discretized solution. */
        double* const a;
        /** A copy of the discretized solution used during the finite
         * difference update. */
        double* const b;
        /** A pointer to the saved solution snapshots. */
        double *z;
        v_visco(int m_,double alpha_);
        ~v_visco();
        void initialize();
        void integrate(double duration,double sfac,bool save,int snaps_);
        void step(double nu);
        void output(const char* filename);
    private:
	inline double eno2(double p0,double p1,double p2,double p3);
        /** The number of snapshots that are saved. */
        int snaps;
};

#endif
